from io import StringIO
import pandas as pd
import requests
import concurrent.futures
import json
from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse
from django.contrib.messages.views import SuccessMessageMixin
from django.urls import reverse, reverse_lazy
from django.db.utils import IntegrityError
from django.core.serializers import serialize

from django.contrib import messages
from django.db.models import ProtectedError, Count, Q, Max, Min
from .models import *
from .forms import *
from django.views.generic import (
    ListView, CreateView, DeleteView, UpdateView, DetailView)
from Bio import SeqIO, Align
from django.db import transaction



def index(request):
    return render(request, "targets/home.html")



# ============== PROJECT ================================

class ProjectListView(ListView):
    template_name_suffix = "_list"
    context_object_name = 'project_list'
    model = Project


class ProjectFormView(SuccessMessageMixin, CreateView):
    model = Project
    template_name_suffix = '_new'
    form_class = ProjectForm
    success_message = "Project was successfully added: %(name)s"

    def get_success_url(self):
        return reverse('targets:project_detail', args=(self.object.id,))


class ProjectDetailView(DetailView):
    model = Project

class ProjectUpdateView(SuccessMessageMixin, UpdateView):
    model = Project
    template_name_suffix = '_update'
    form_class = ProjectForm
    success_message = "Project was successfully updated:  %(name)s"
    
    def get_success_url(self):
        return reverse('targets:project_detail', args=(self.object.id,))
    
    def get_initial(self):
        initial = super().get_initial()
        # Retrieve the object being updated
        obj = self.get_object()
        # Populate the initial data with values from the object
        initial['name'] = obj.name
        initial['description'] = obj.description
        initial['notes'] = obj.notes
        return initial


class ProjectDeleteView(DeleteView):
    model = Project
    success_url = reverse_lazy('targets:project_list', args=())
    
    def post(self, request, *args, **kwargs):
        try:
            return self.delete(request, *args, **kwargs)
        except ProtectedError:
            return render(request, "targets/protected_error.html")




# ============== ORGANISM ================================


def select_organism_list(request, pk):
    if request.method == 'POST':
        organism_resp = request.POST
        organism, created = Organism.objects.update_or_create(
                name=organism_resp.get('species'),
                defaults={
                    'db_isolates': organism_resp.get('isolates'),
                    'db_seqdef':organism_resp.get('seqdef')}
            )
        project = Project.objects.get(id=pk)
        # remove genomes and loci if organism is changed
        if (project.organism) and (project.organism.id != organism.id):
            project.loci.clear()
            project.genomes.clear()
        # Save any updates to the database fields
        project.organism = organism
        project.save()
        if created:
            messages.success(
                request, f"Successfully created {organism.name} created and added to project {project.name}.")
        else:
            messages.success(
                request, f"Succesfully updated {organism.name} and added to project {project.name}.")
        return redirect('targets:project_detail', pk=pk)
        
    return render(request, "targets/select_organism.html",  {'pk': pk})


class OrganismListView(ListView):
    template_name_suffix = "_list"
    context_object_name = 'organism_list'
    model = Organism


class OrganismDetailView(DetailView):
    model = Organism


def get_organism_json(request):
    r = requests.get('https://rest.pubmlst.org/db', params=request.GET)
    res = json.loads(r.text)
    data = []
    idx = 0
    org_dict = {}
    for genus in res:
        for db in genus['databases']:
            if db['name'].split("_")[-1] == 'isolates':
                species = db['description'].replace(' isolates', '')
                if species in org_dict:
                    org_dict[species]['isolates'] = db['name']
                else:
                    org_dict[species] = {'isolates': db['name']}
            if db['name'].split("_")[-1] == "seqdef":
                species = db['description'].replace(' sequence/profile definitions', '')
                if species in org_dict:
                    org_dict[species]['seqdef'] = db['name']
                else:
                    org_dict[species] = {'seqdef': db['name']}
    for org, dbs in org_dict.items():
        if len(dbs) > 1: # don't allow species with no seqdef
            data.append({
                'id': idx, 
                'species': org,
                'isolates': dbs['isolates'],
                'seqdef': dbs['seqdef']})
            idx += 1
    return JsonResponse(data, safe=False)

# class OrganismUpdateView(UpdateView):
#     model = Panel
#     template_name_suffix = '_update'
#     form_class = PanelForm
#     success_message = "Panel was successfully updated:  %(name)s"
    
#     def get_success_url(self):
#         return reverse('targets:project_detail', args=(self.object.id,))


# ============== Data ==================================


def load_project_data(request, pk):
    if request.method == 'POST':
        form = DataUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file_name = request.FILES['file']
            project = Project.objects.get(pk=pk)
            # Set all columns and row to inactive
            project.deactivate_all_genomes()
            project.deactivate_all_loci()
            project.deactivate_all_metadata_categories()

            organism = project.organism
            
            df = pd.read_excel(file_name, index_col='id')
            print(df.head())
            # pull existing genome and loci data from database
            # to filter out any values from the table that may
            # not belong.
            db_genomes = organism.get_genomes_list(request)
            db_loci = organism.get_loci_list(request)
    
            loci_cols = df.columns.intersection(db_loci)
            genome_idx = df.index.intersection(db_genomes)

            threshold_percentage = 75

            # Calculate the minimum number of non-NaN values required based on the threshold percentage
            threshold_count = int((100 - threshold_percentage) / 100 * len(loci_cols))

            # Drop rows with more than the specified threshold count of NaN values
            df = df.dropna(thresh=threshold_count)

            # Consider all other columns to be metadata
            metadata_cols = df.columns.difference(db_loci)
            # keep only rows (genomes) that are found in the database
            meta_df = df.loc[genome_idx, metadata_cols].to_dict()
            alleles_df = df.loc[genome_idx, loci_cols].to_dict()
            # TODO: add warning if alleles/genomes are removed
            # add loci to database
            loci_dict = {}
            print("creating LOCI")
            for locus in loci_cols:
                locus_inst, created = Locus.objects.get_or_create(
                    name=locus,
                    project=project,
                )
                if not created:
                    # Created are active by default but if
                    # locus already exists move to active
                    locus_inst.active = True
                    locus_inst.save()
                loci_dict[locus] = locus_inst
            print("creating GENOMES")

            # add genomes to database
            genome_dict = {}
            for genome in genome_idx:
                genome_inst, created = Genome.objects.get_or_create(
                    name=genome,
                    project=project,
                )
                if not created:
                    # Created are active by default but if
                    # genome already exists move to active
                    genome_inst.active = True
                    genome_inst.save()
                genome_dict[genome] = genome_inst
            print("creating metacat")

            # add metadata cateories to database
            meta_dict = {}
            for meta in metadata_cols:
                meta_inst, created = MetadataCategory.objects.get_or_create(
                    name=meta,
                    project=project,
                )
                if not created:
                    # Created are active by default but if
                    # metadata cat already exists move to active
                    meta_inst.active = True
                    meta_inst.save()
                meta_dict[meta] = meta_inst

            # add alleles
            print("creating alleles")

            for locus, values in alleles_df.items():
                locus_inst = loci_dict[locus]
                for genome, allele in values.items():
                    if str(allele) != 'nan': # don't add alleles with missing data
                        Allele.objects.update_or_create(
                            genome=genome_dict[genome],
                            locus=locus_inst,
                            project=project,
                            defaults={'value': allele}
                            )
            print("creating METADATA")

            # add metadata
            for metacat, values in meta_df.items():
                metacat_inst = meta_dict[metacat]
                for genome, value in values.items():
                    if str(value) != 'nan': # don't add metadata with missing values
                        Metadata.objects.update_or_create(
                            genome=genome_dict[genome],
                            category=metacat_inst,
                            project=project,
                            defaults={'value': value}
                            )
                       

            return redirect('targets:project_detail', pk=pk)
    else:
        form = DataUploadForm()

    return render(request, 'targets/data_upload.html', {'form': form})



# ============== GENOME ==================================
class GenomeListView(ListView):
    template_name = "targets/genome_list.html"
    model = Genome
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_active_genomes()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_genomes_json(request, pk, active=None):
    project = Project.objects.get(id=pk)
    if active is None:
        genomes = project.get_genomes()
    else:
        genomes = project.get_active_genomes()

    all_rows = []

    for genome in genomes:
        all_rows.append({
            'id': genome.id,
            'name': genome.name,
            'alleles': genome.count_alleles(),
            'metadata': genome.count_metadata(),
            'active': genome.active
        })


    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)




def select_genomes_list(request, pk):
    if request.method == 'POST':
        genome_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deactivate_all_genomes()
        activate = Genome.objects.filter(project=project, pk__in=genome_resp)
        for genome in activate:
            genome.active = True
            genome.save()
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_genomes.html", {'pk': pk})




# ============== LOCI ==================================
class LocusListView(ListView):
    template_name = "targets/locus_list.html"
    model = Locus
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_active_loci()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_loci_json(request, pk, active=None):
    project = Project.objects.get(id=pk)
    if active is None:
        loci = project.get_loci()
    else:
        loci = project.get_active_loci()

    all_rows = []
    for locus in loci:
        all_rows.append({
            'id': locus.id,
            'name': locus.name,
            'alleles': locus.count_alleles,
            'active': locus.active
        })

    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)


def select_loci_list(request, pk):
    if request.method == 'POST':
        loci_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deactivate_all_loci()
        activate = Locus.objects.filter(project=project, pk__in=loci_resp)
        for locus in activate:
            locus.active = True
            locus.save()
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_loci.html", {'pk': pk})


# ============== METADATA ==================================
class MetadataCategoryListView(ListView):
    template_name = "targets/metadatacategory_list.html"
    model = MetadataCategory
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_active_metadata()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_metacat_json(request, pk, active=None):
    project = Project.objects.get(id=pk)
    if active is None:
        meta = project.get_metadata()
    else:
        meta = project.get_active_metadata()

    all_rows = []
    for dat in meta:
        all_rows.append({
            'id': dat.id,
            'name': dat.name,
            'values': dat.count_values,
            'active': dat.active
        })

    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)


def select_metacat_list(request, pk):
    if request.method == 'POST':
        meta_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deactivate_all_metadata_categories()
        activate = MetadataCategory.objects.filter(project=project, pk__in=meta_resp)
        for meta in activate:
            meta.active = True
            meta.save()
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_metadatacategory.html", {'pk': pk})

# ============== DATA ==================================


def data_list_view(request, pk):
    project = Project.objects.get(id=pk)
    active_genomes = project.get_active_genomes()
    active_loci = project.get_active_loci()
    active_metadata = project.get_active_metadata()
    alleles = Allele.objects.filter(
        project=project,
        genome__in=active_genomes.values_list('id', flat=True),
        locus__in=active_loci.values_list('id', flat=True)
        ).values('genome', 'locus', 'value')
    meta = Metadata.objects.filter(
        project=project,
        genome__in=active_genomes.values_list('id', flat=True),
        category__in=active_metadata.values_list('id', flat=True)
        ).values('genome', 'category', 'value')
    allele_df = pd.DataFrame.from_dict(alleles)
    allele_df = allele_df.pivot(index="genome", columns='locus', values='value')
    idx_dict = {gen.id: gen.name for gen in active_genomes}
    allele_df.index = [idx_dict.get(i) for i in allele_df.index]
    col_dict = {loc.id: loc.name for loc in active_loci}
    allele_df.columns = [col_dict.get(i) for i in allele_df.columns]

    meta_df = pd.DataFrame.from_dict(meta)
    meta_df = meta_df.pivot(index="genome", columns='category', values='value')
    idx_dict = {gen.id: gen.name for gen in active_genomes}
    meta_df.index = [idx_dict.get(i) for i in meta_df.index]
    col_dict = {loc.id: loc.name for loc in active_metadata}
    meta_df.columns = [col_dict.get(i) for i in meta_df.columns]

    df = meta_df.join(allele_df)
    df = df.fillna('')
    print(df.head())
    html = df.to_html(
        table_id="datatable",
        classes=['display', 'table', 'table-hover'],
        justify="left", index=True)
    print(html[:100])

    return render(
                request, 'targets/data_table.html',
                {'table': html, 'project': project}
                )

    

def data_delete_view(request, pk):
    project = Project.objects.get(pk=pk)
    if request.method == 'POST':
        project.get_genomes().delete()
        project.get_loci().delete()
        project.get_metadata().delete()
        # Handle the confirmation of the delete action
        return redirect('targets:project_list')  # Redirect to a success page after deletion

    # Render the confirmation template for the delete action
    return render(request, 'targets/data_confirm_delete.html', {'project': project})

# ============== LOCI ==================================

# class LocusListView(ListView):
#     template_name = "targets/locus_list.html"
#     model = Locus
#     def get_queryset(self, **kwargs):
#         new_context = Project.objects.get(pk=self.kwargs['pk']).loci.all()
#         return new_context
    
#     def get_context_data(self, **kwargs):
#         context = super().get_context_data(**kwargs)
#         project = Project.objects.get(pk=self.kwargs['pk'])
#         context['project'] = project
#         return context
    


# def get_loci_json(request, pk):
#     project = Project.objects.get(id=pk)
#     r = requests.get(
#         f'https://rest.pubmlst.org/db/{project.organism.db_isolates}/loci?return_all=1', params=request.GET)
#     split = str.split
#     res = json.loads(r.text)
#     data = [
#         {'id': i, 'description': v} for i, v in enumerate(
#         map(lambda x: split(x, "/")[-1], res['loci']))]

#     return JsonResponse(data, safe=False)


# def select_loci_list(request, pk):
#     if request.method == 'POST':
#         loci_resp = request.POST
#         created_count = 0
#         updated_count = 0
#         project = Project.objects.get(id=pk)
#         project.loci.clear()
#         add_loci = []
#         for locus in loci_resp.getlist('loci'):
#             locus, created = Locus.objects.get_or_create(
#                     name=locus,
#                     organism=project.organism,
#                 )
#             add_loci.append(locus)
#             if created:
#                 created_count += 1
#             else:
#                 updated_count += 1
#         project.loci.add(*add_loci)
#         project.save()
#         msg = f"Successfully added {created_count} loci project {project.name}. " if created_count else ""
#         msg += f"Successfully updated {updated_count} loci for project {project.name}" if updated_count else ""
#         messages.success(request, msg)
#         return redirect('targets:project_detail', pk=pk) # TODO: change redirect to loci list when created
        
#     return render(request, "targets/select_loci.html", {'pk': pk})
#     # r = requests.get('https://rest.pubmlst.org/db', params=request.GET)


# # ============== PANELS ================================


# class PanelListView(ListView):
#     template_name_suffix = "_list"
#     context_object_name = 'panel_list'
#     model = Panel

# class PanelUpdateView(SuccessMessageMixin, UpdateView):
#     model = Panel
#     template_name_suffix = '_update'
#     form_class = PanelForm
#     success_message = "Panel was successfully updated:  %(name)s"
    
#     def get_success_url(self):
#         return reverse('targets:project_detail', args=(self.object.id,))


# class PanelDetailView(DetailView):
#     model = Panel


# class PanelDeleteView(DeleteView):
#     model = Panel
#     success_url = reverse_lazy('targets:panel_list', args=())
    
#     def post(self, request, *args, **kwargs):
#         try:
#             return self.delete(request, *args, **kwargs)
#         except ProtectedError:
#             return render(request, "targets/protected_error.html")






# def get_alignment(ref, seq):
#     aligner = Align.PairwiseAligner()
#     aligner.mode = 'global'
#     alignments = aligner.align(ref, seq)










# def get_loci_json(request, organism):
#     org = Organism.objects.get(id=organism)
#     r = requests.get(
#         'https://rest.pubmlst.org/db/{database}/loci?return_all=1'.format(
#             database=org.isolates), params=request.GET)
#     split = str.split
#     res = json.loads(r.text)
#     data = [
#         {'id': i, 'description': v} for i, v in enumerate(
#         map(lambda x: split(x, "/")[-1], res['loci']))]

#     return JsonResponse(data, safe=False)



# class OrganismLociListView(ListView):
#     template_name = "targets/organism_loci_list.html"
#     model = Locus

#     def get_queryset(self, **kwargs):
#         new_context = Locus.objects.filter(
#             organism=self.kwargs['pk'],
#         )
#         return new_context
    
