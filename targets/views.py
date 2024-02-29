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
from django.core.files.storage import default_storage
from django.core.files import File
from django.core.files.base import ContentFile


from .utils import *


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
                request, f"Successfully created {organism.name} and added to project {project.name}.")
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



# ============== GENOME ==================================
class GenomeListView(ListView):
    template_name = "targets/genome_list.html"
    model = Genome
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_selected_genomes()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_genomes_json(request, pk, selected=False):
    project = Project.objects.get(id=pk)
    if selected:
        genomes = project.get_selected_genomes()
    else:
        genomes = project.get_active_genomes()

    all_rows = []

    for genome in genomes:
        all_rows.append({
            'id': genome.id,
            'name': genome.name,
            'alleles': genome.n_alleles,
            'metadata': genome.n_metadata,
            'active': genome.active,
            'selected': genome.selected
        })


    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)




def select_genomes_list(request, pk):
    if request.method == 'POST':
        genome_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deselect_all_genomes()
        Genome.objects.filter(
            project=project, pk__in=genome_resp).update(selected=True)
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_genomes.html", {'pk': pk})




# ============== LOCI ==================================
class LocusListView(ListView):
    template_name = "targets/locus_list.html"
    model = Locus
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_selected_loci()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_loci_json(request, pk, selected=False):
    project = Project.objects.get(id=pk)
    if selected:
        loci = project.get_selected_loci()
    else:
        loci = project.get_active_loci()

    all_rows = []
    for locus in loci:
        all_rows.append({
            'id': locus.id,
            'name': locus.name,
            'alleles': locus.n_alleles,
            'active': locus.active,
            'selected': locus.selected
        })

    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)


def select_loci_list(request, pk):
    if request.method == 'POST':
        loci_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deselect_all_loci()
        Locus.objects.filter(project=project, pk__in=loci_resp).update(selected=True)
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_loci.html", {'pk': pk})


# ============== METADATA ==================================
class MetadataCategoryListView(ListView):
    template_name = "targets/metadatacategory_list.html"
    model = MetadataCategory
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).get_selected_metadata()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_metacat_json(request, pk, selected=False):
    project = Project.objects.get(id=pk)
    if selected:
        meta = project.get_selected_metadata()
    else:
        meta = project.get_active_metadata()

    all_rows = []
    for dat in meta:
        all_rows.append({
            'id': dat.id,
            'name': dat.name,
            'values': dat.n_values,
            'genomes': dat.n_genomes,
            'active': dat.active,
            'selected': dat.selected
        })

    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)


def select_metacat_list(request, pk):
    if request.method == 'POST':
        meta_resp = map(int, request.POST.getlist('id'))
        project = Project.objects.get(pk=pk)
        project.deselect_all_metadata()
        MetadataCategory.objects.filter(project=project, pk__in=meta_resp).update(selected=True)
        return redirect('targets:project_detail', pk=pk)

    return render(request, "targets/select_metadatacategory.html", {'pk': pk})



# def load_project_data(request, pk):
#     if request.method == 'POST':
#         form = DataUploadForm(request.POST, request.FILES)
#         if form.is_valid():
#             file_name = request.FILES['file']
#             project = Project.objects.get(pk=pk)
#             # Remove the existing file (if any)
#             if project.data:
#                 project.data.delete()
#             if project.imputed_data:
#                 project.imputed_data.delete()
#             project.data = file_name
#             project.save()
#             # Read data from file
#             try:
#                 df = pd.read_excel(project.data, index_col='id')
#             except pd.errors.ExcelError:
#                 project.data = None
#                 project.save()
#                 messages.error(
#                     request,
#                     'Failed to load data, please check formatting and retry.',
#                     extra_tags="danger")
#                 return redirect('targets:project_detail', pk=pk)
#             try:
#                 organism = project.organism
#                 # pull existing loci data from database
#                 # to filter out any values from the table that may
#                 # not belong. 
#                 db_loci = organism.get_loci_list(request)
#                 loci_cols = df.columns.intersection(db_loci)

#                 if len(loci_cols) == 0:
#                     raise IOError("Loci not found for selected organism")
#             except AttributeError:
#                 project.data = None
#                 project.save()
#                 messages.error(
#                     request, 'No Organism! Select organism before uploading data.',
#                     extra_tags="danger")
#                 return redirect('targets:project_detail', pk=pk)
#             except IOError:
#                 project.data=None
#                 project.save()
#                 messages.error(
#                     request,
#                     "None of the loci are found for organism. Please check that data matches selected organism.",
#                     extra_tags="danger")
#                 return redirect('targets:project_detail', pk=pk)
            
#             msg = f"Data successfully added to project. " # add messages 



#             # Consider all other columns to be metadata
#             genome_idx = df.index
#             metadata_cols = df.columns.difference(db_loci)
#             meta_df = df.loc[genome_idx, metadata_cols]
#             alleles_df = df.loc[genome_idx, loci_cols]

#             # Drop Genomes (rows) with more than the specified threshold of NaN values
            
#             alleles_df = drop_nans(alleles_df, 0.8)
#             # record dropped genomes
#             n_genomes_dropped = len(genome_idx.difference(alleles_df.index))
#             genome_idx = alleles_df.index # update to reflect dropped rows
#             alleles_df = alleles_df.loc[genome_idx, loci_cols]

#             # Drop Loci (cols) with more than the specified threshold of NaN values
#             alleles_df = drop_nans(alleles_df, 0.8, axis=1)
#             # record dropped alleles
#             n_loci_dropped = len(loci_cols.difference(alleles_df.columns))
#             # update loci columns to reflect dropped 
#             loci_cols = alleles_df.columns
#             # alleles_df = df.loc[genome_idx, loci_cols]

#             print("LOCI Dropped", n_loci_dropped)
#             print("Genomes Dropped", n_genomes_dropped)
#             if n_loci_dropped:
#                 msg += f"{n_loci_dropped} loci were dropped due to missing data. " 
#             if n_genomes_dropped:
#                 msg += f"{n_genomes_dropped} genomes were dropped due to missing data. "
#             # add LOCI to database
#             print("UPDATING Loci")

#             # ids of existing loci that should be activated (in uploaded table)
#             Locus.objects.filter(
#                 project=project,
#                 name__in=loci_cols).update(active=True)
#             # ids of existing loci that should be deactivated (not in uploaded table)
#             Locus.objects.filter(project=project).exclude(
#                 name__in=loci_cols).update(active=False)
#             # New loci that are not in database to create
#             new_loci = list(set(loci_cols) - set(Locus.objects.filter(
#                 project=project,
#                 name__in=loci_cols).values_list('name', flat=True)))
#             # bulk create new loci
#             Locus.objects.bulk_create(
#                 [Locus(project=project, name=name) for name in new_loci])
#             # count number of unique alleles and save to instance
#             all_loci = Locus.objects.filter(project=project, active=True)
#             for locus in all_loci:
#                 locus.alleles = alleles_df[locus.name].dropna().nunique()
#                 locus.save()

#             # add GENOMES to database
#             print("UPDATING genomes")
#             # ids of existing genomes that should be activated (in uploaded table)
#             Genome.objects.filter(
#                 project=project,
#                 name__in=genome_idx).update(active=True)
#             # ids of existing genomes that should be deactivated (not in uploaded table)
#             Genome.objects.filter(project=project).exclude(
#                 name__in=genome_idx).update(active=False)
#             # New Genomes that are not in database to create
#             new_genomes = list(set(genome_idx) - set(Genome.objects.filter(
#                 project=project,
#                 name__in=genome_idx).values_list('name', flat=True)))
#             # bulk create new genomes
#             Genome.objects.bulk_create(
#                 [Genome(project=project, name=name) for name in new_genomes])
            
#             # count number of unique alleles and save to instance
#             all_genomes = Genome.objects.filter(project=project, active=True)
#             for genome in all_genomes:
#                 genome.alleles = alleles_df.loc[genome.name].dropna().size
#                 genome.metadata = meta_df.loc[genome.name].dropna().size
#                 genome.save()
            

#             # add METADATA CATEGORIES to database
#             print("UPDATING metacategories")
#             # ids of existing metacats that should be activated (in uploaded table)
#             MetadataCategory.objects.filter(
#                 project=project,
#                 name__in=metadata_cols).update(active=True)
#             # ids of existing metacats that should be deactivated (not in uploaded table)
#             MetadataCategory.objects.filter(project=project).exclude(
#                 name__in=metadata_cols).update(active=False)
#             # New metacats that are not in database to create
#             new_metacat = list(set(metadata_cols) - set(MetadataCategory.objects.filter(
#                 project=project,
#                 name__in=metadata_cols).values_list('name', flat=True)))
#             # bulk create new metacats
#             MetadataCategory.objects.bulk_create(
#                 [MetadataCategory(project=project, name=name) for name in new_metacat])

#             # count number of genomes with data and unique values and save to instance
#             all_meta = MetadataCategory.objects.filter(project=project, active=True)
#             for meta in all_meta:
#                 meta.values = meta_df[meta.name].dropna().nunique()
#                 meta.genomes = meta_df[meta.name].dropna().size
#                 meta.save()

#             # impute data
#             print("IMPUTING DATA")
#             imputed_data = impute_missing_alleles(project.get_allele_table())
#             print(imputed_data)
#             # Convert the DataFrame to CSV format as a string
#             csv_content = imputed_data.to_csv(index=False, sep="\t")

#             # Create a ContentFile from the CSV content
#             content_file = ContentFile(csv_content)

#             # Create a File instance from the ContentFile
#             file_instance = File(content_file)

#             # Create an instance of YourModel and associate the File instance with its file_field
#             project.imputed_data.save(f'imputed_data_{project.id}.tsv', file_instance)
#             project.save()
#             messages.success(request, msg)
#             return redirect('targets:project_detail', pk=pk)
#     else:
#         form = DataUploadForm()

#     return render(request, 'targets/data_upload.html', {'form': form})

# def data_list_view(request, pk):
#     project = Project.objects.get(id=pk)
#     allele_table = project.get_allele_table()
#     metadata_table = project.get_metadata_table()
#     df = metadata_table.join(allele_table)
#     df.index.name = "Genome ID"
#     html = df.reset_index().to_html(
#         table_id="datatable",
#         classes=['display', 'table', 'table-hover'],
#         justify="left", index=False)
#     metadata_idx = [i for i in range(1, metadata_table.shape[1] + 1)]
#     allele_idx = [i for i in range(metadata_table.shape[1] + 1, df.shape[1] + 1)]

#     return render(
#                 request, 'targets/data_table.html',
#                 {
#                     'table': html,
#                     'project': project,
#                     'meta_idx': metadata_idx,
#                     'allele_idx': allele_idx
#                 }
#                 )



def load_project_data(request, pk):
    if request.method == 'POST':
        form = DataUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file_name = request.FILES['file']
            project = Project.objects.get(pk=pk)
            
            # Read data from file
            try:
                df = pd.read_excel(file_name, index_col=0, dtype=str, na_values='',
                                   keep_default_na=False).fillna('')
                df.index.name="#Genome"
            except pd.errors.ExcelError:
                try:
                    df = pd.read_csv(file_name, index_col=0, dtype=str, na_values='',
                                   keep_default_na=False).fillna('')
                    df.index.name="#Genome"
                except Exception:
                    messages.error(
                        request,
                        'Failed to load data, please check formatting and retry.',
                        extra_tags="danger")
                    return redirect('targets:project_detail', pk=pk)
            try:
                organism = project.organism
                # pull existing loci data from database
                # to filter out any values from the table that may
                # not belong. 
                db_loci = organism.get_loci_list(request)
                loci_cols = df.columns.intersection(db_loci)

                if len(loci_cols) == 0:
                    raise IOError("Loci not found for selected organism")
            except AttributeError:
                messages.error(
                    request, 'No Organism! Select organism before uploading data.',
                    extra_tags="danger")
                return redirect('targets:project_detail', pk=pk)
            except IOError:
                messages.error(
                    request,
                    "None of the loci are found for organism. Please check that data matches selected organism.",
                    extra_tags="danger")
                return redirect('targets:project_detail', pk=pk)
            
            msg = f"Data successfully added to project. " # add messages 



            # Consider all other columns to be metadata
            genome_idx = df.index
            metadata_cols = df.columns.difference(db_loci)
            meta_df = df.loc[genome_idx, metadata_cols]
            alleles_df = df.loc[genome_idx, loci_cols]

            # Drop Genomes (rows) with more than the specified threshold of NaN values
            
            alleles_df = drop_nans(alleles_df, 0.8)
            # record dropped genomes
            n_genomes_dropped = len(genome_idx.difference(alleles_df.index))
            genome_idx = alleles_df.index # update to reflect dropped rows
            alleles_df = alleles_df.loc[genome_idx, loci_cols]

            # Drop Loci (cols) with more than the specified threshold of NaN values
            alleles_df = drop_nans(alleles_df, 0.8, axis=1)
            # record dropped alleles
            n_loci_dropped = len(loci_cols.difference(alleles_df.columns))
            # update loci columns to reflect dropped 
            loci_cols = alleles_df.columns
            alleles_df = df.loc[genome_idx, loci_cols]

            print("LOCI Dropped", n_loci_dropped)
            print("Genomes Dropped", n_genomes_dropped)
            if n_loci_dropped:
                msg += f"{n_loci_dropped} loci were dropped due to missing data. " 
            if n_genomes_dropped:
                msg += f"{n_genomes_dropped} genomes were dropped due to missing data. "
            # add LOCI to database
            print("UPDATING Loci")

            # ids of existing loci that should be activated (in uploaded table)
            Locus.objects.filter(
                project=project,
                name__in=loci_cols).update(active=True, selected=True)
            # ids of existing loci that should be deactivated (not in uploaded table)
            Locus.objects.filter(project=project).exclude(
                name__in=loci_cols).update(active=False, selected=False)
            # New loci that are not in database to create
            new_loci = list(set(loci_cols) - set(Locus.objects.filter(
                project=project,
                name__in=loci_cols).values_list('name', flat=True)))
            # bulk create new loci
            Locus.objects.bulk_create(
                [Locus(project=project, name=name) for name in new_loci])
            # count number of unique alleles and save to instance
            all_loci = Locus.objects.filter(project=project, active=True, selected=True)
            for locus in all_loci:
                locus.n_alleles = alleles_df[locus.name].dropna().nunique()
                locus.save()

            # add GENOMES to database
            print("UPDATING genomes")
            # ids of existing genomes that should be activated (in uploaded table)
            Genome.objects.filter(
                project=project,
                name__in=genome_idx).update(active=True, selected=True)
            # ids of existing genomes that should be deactivated (not in uploaded table)
            Genome.objects.filter(project=project).exclude(
                name__in=genome_idx).update(active=False, selected=False)
            # New Genomes that are not in database to create
            new_genomes = list(set(genome_idx) - set(Genome.objects.filter(
                project=project,
                name__in=genome_idx).values_list('name', flat=True)))
            # bulk create new genomes
            Genome.objects.bulk_create(
                [Genome(project=project, name=name) for name in new_genomes])
            
            # count number of unique alleles and save to instance
            all_genomes = Genome.objects.filter(project=project, active=True, selected=True)
            for genome in all_genomes:
                genome.n_alleles = alleles_df.loc[genome.name].dropna().size
                genome.n_metadata = meta_df.loc[genome.name].dropna().size
                genome.save()
            

            # add METADATA CATEGORIES to database
            print("UPDATING metacategories")
            # ids of existing metacats that should be activated (in uploaded table)
            MetadataCategory.objects.filter(
                project=project,
                name__in=metadata_cols).update(active=True, selected=True)
            # ids of existing metacats that should be deactivated (not in uploaded table)
            MetadataCategory.objects.filter(project=project).exclude(
                name__in=metadata_cols).update(active=False, selected=False)
            # New metacats that are not in database to create
            new_metacat = list(set(metadata_cols) - set(MetadataCategory.objects.filter(
                project=project,
                name__in=metadata_cols).values_list('name', flat=True)))
            # bulk create new metacats
            MetadataCategory.objects.bulk_create(
                [MetadataCategory(project=project, name=name) for name in new_metacat])

            # count number of genomes with data and unique values and save to instance
            all_meta = MetadataCategory.objects.filter(project=project, active=True, selected=True)
            for meta in all_meta:
                meta.n_values = meta_df[meta.name].dropna().nunique()
                meta.n_genomes = meta_df[meta.name].dropna().size
                meta.save()

            # impute data
            print("IMPUTING DATA")
            imputed_data = impute_missing_alleles(alleles_df).set_index('Genome')
            alleles_df.index.name = 'genome'
            alleles_df = pd.melt(
                alleles_df, ignore_index=False,
                var_name = "locus", value_name="allele")
            alleles_df = alleles_df.reset_index().set_index(['genome', 'locus'])
            imputed_data.index.name = 'genome'
            imputed_data.columns = ['locus', 'imputed']

            imputed_data = imputed_data.reset_index().set_index(['genome', 'locus'])
            alleles_df = alleles_df.join(imputed_data).reset_index()
            alleles_df['project'] = project
            # replace locus names with ids
            locus_map = {k: Locus.objects.get(name=k) for k in alleles_df['locus'].unique()}
            alleles_df['locus'] = alleles_df['locus'].replace(locus_map)
            genome_map = {k: Genome.objects.get(name=k) for k in alleles_df['genome'].unique()}
            alleles_df['genome'] = alleles_df['genome'].replace(genome_map)
            alleles_df = alleles_df.set_index(['project', 'genome', 'locus'])
            # Adding Alleles
            print("ADDING ALLELES")
            Allele.objects.bulk_create(
                [Allele(
                    project=dat[0],
                    genome=dat[1],
                    locus=dat[2],
                    allele=dat[3],
                    imputed=dat[4]
                    ) for dat in alleles_df.to_records()], 
                update_conflicts=True,
                unique_fields=['project', 'genome', 'locus'],
                update_fields=['allele', 'imputed'],
            )
            # select metadata using current genome idx
            print("ADDING METADATA")
            meta_df = meta_df.loc[genome_idx, metadata_cols]
            meta_df.index.name = 'genome'
            meta_df = pd.melt(
                meta_df, ignore_index=False,
                var_name='category', value_name='value')
            meta_df['project'] = project
            meta_df = meta_df.reset_index()
            cat_map = {k: MetadataCategory.objects.get(name=k) for k in meta_df['category'].unique()}
            meta_df['category'] = meta_df['category'].replace(cat_map)
            genome_map = {k: Genome.objects.get(name=k) for k in meta_df['genome'].unique()}
            meta_df['genome'] = meta_df['genome'].replace(genome_map)
            meta_df = meta_df.set_index(['project', 'genome', 'category'])
            Metadata.objects.bulk_create(
                [Metadata(
                    project=dat[0],
                    genome=dat[1],
                    category=dat[2],
                    value=dat[3]
                ) for dat in meta_df.to_records()], 
                update_conflicts=True,
                unique_fields=['project', 'genome', 'category'],
                update_fields=['value'],
            )


            messages.success(request, msg)
            return redirect('targets:project_detail', pk=pk)
    else:
        form = DataUploadForm()

    return render(request, 'targets/data_upload.html', {'form': form})

def data_list_view(request, pk):
    project = Project.objects.get(id=pk)
    allele_table = project.get_allele_table()
    metadata_table = project.get_metadata_table()
    print(allele_table)
    
    df = metadata_table.join(allele_table)
    df.index.name = "Genome ID"
    html = df.reset_index().to_html(
        table_id="datatable",
        classes=['display', 'table', 'table-hover'],
        justify="left", index=False)
    metadata_idx = [i for i in range(1, metadata_table.shape[1] + 1)]
    allele_idx = [i for i in range(metadata_table.shape[1] + 1, df.shape[1] + 1)]

    return render(
                request, 'targets/data_table.html',
                {
                    'table': html,
                    'project': project,
                    'meta_idx': metadata_idx,
                    'allele_idx': allele_idx
                }
                )


def data_delete_view(request, pk):
    project = Project.objects.get(pk=pk)
    if request.method == 'POST':

        project.deactivate_all_genomes()
        project.deactivate_all_loci()
        project.deactivate_all_metadata_categories()
        # Delete all targetsets associated with the data
        # TargetCollection.objects.filter(project=project).delete()

    
        # Handle the confirmation of the delete action
        return redirect('targets:project_list')  # Redirect to a success page after deletion

    # Render the confirmation template for the delete action
    return render(request, 'targets/data_confirm_delete.html', {'project': project})

# ============== TARGET COLLECTION ==================================



class TargetCollectionFormView(CreateView):
    model = TargetCollection
    template_name_suffix = '_new'
    form_class = TargetCollectionForm
    success_message = "Target collection was successfully added: %(name)s"

    def get_success_url(self):
        return reverse('targets:target_collection_detail', args=(self.object.id,))
    
    def get_initial(self):
        initial = super().get_initial()
        print(self.kwargs.get('pk'))
        initial['project'] = self.kwargs.get('pk')
        return initial


class TargetCollectionDetailView(DetailView):
    model = TargetCollection
 

class TargetCollectionUpdateView(SuccessMessageMixin, UpdateView):
    model = TargetCollection
    template_name_suffix = '_update'
    form_class = TargetCollectionUpdateForm
    success_message = "Target collection was successfully updated:  %(name)s"

    
    def get_success_url(self):
        return reverse('targets:target_collection_detail', args=(self.object.id,))
    
    def get_initial(self):
        initial = super().get_initial()
        # Retrieve the object being updated
        obj = self.get_object()
        # Populate the initial data with values from the object
        initial['project'] = obj.project.id
        initial['name'] = obj.name
        initial['description'] = obj.description
        return initial



class TargetCollectionListView(ListView):
    template_name_suffix = "_list"
    model = TargetCollection

    def get_queryset(self):
        # Access kwargs using self.kwargs
        project_pk = self.kwargs.get('pk', None)
        print(project_pk)
        # Filter by project if pk is passed
        if project_pk:
            collection = TargetCollection.objects.filter(project=project_pk)
        else:
            collection = TargetCollection.objects.all()
        return collection
    
    def get_context_data(self):
        """
        Add context for add targetcollection button,
        sets project pk as initial if in project based list view.
        """
        context = super().get_context_data(**self.kwargs)
        if self.kwargs.get('pk', None):
            context['pk'] = self.kwargs.get('pk', None)
        return context
            


class TargetCollectionDeleteView(DeleteView):
    model = TargetCollection

    def get_success_url(self):
        return reverse('targets:target_collection_list')
    
    def post(self, request, *args, **kwargs):
        try:
            return self.delete(request, *args, **kwargs)
        except ProtectedError:
            return render(request, "targets/protected_error.html")
        

class TargetCollectionLocusListView(DetailView):
    template_name = "targets/targets_list.html"
    model = TargetCollection

    

def get_targets_json(request, pk, selected=False):
    """
    Return either all active loci for projects with 
    selected targets marked [selected=False] or
    the targets that have been selected for the 
    target collection [selected=True]
    """
    target_col = TargetCollection.objects.get(id=pk)
    project = target_col.project
    all_rows = []
    selected_targets = target_col.loci.all()
    for locus in selected_targets:
        all_rows.append({
            'id': locus.id,
            'name': locus.name,
            'alleles': locus.n_alleles,
            'active': True
    })
    # If selected is not set, return all active loci.
    if not selected:
        all_active = project.get_selected_loci().difference(selected_targets)
        for locus in all_active:
            all_rows.append({
                'id': locus.id,
                'name': locus.name,
                'alleles': locus.n_alleles,
                'active': False
            })
    data = {'data': all_rows}
    # Return JSON response
    return JsonResponse(data, safe=False)


def select_targets_manual(request, pk):
    target_col = TargetCollection.objects.get(pk=pk)
    if request.method == 'POST':
        loci_resp = map(int, request.POST.getlist('id'))
        target_col = TargetCollection.objects.get(pk=pk)
        target_col.loci.clear()
        target_col.loci.add(*loci_resp)
        return redirect('targets:target_collection_detail', pk=pk)

    return render(
        request, "targets/select_targets.html",
        {'targetcollection': target_col})



class TargetCollectionAddTargets(SuccessMessageMixin, UpdateView):
    model = TargetCollection
    template_name_suffix = '_add_targets'
    form_class = AddTargets
    success_message = "Targets were successfully added"

    def get_success_url(self):
        return reverse('targets:target_collection_detail', args=(self.object.id,))
    
    def form_valid(self, form):
        
        response = super().form_valid(form)
        target_collection = TargetCollection.objects.get(id=self.kwargs.get('pk'))
        project = target_collection.project
        # alleles_df = project.get_allele_table()
        imputed_df = project.get_imputed_data() 
        print(imputed_df.shape)
        selected_loci = target_collection.loci.all()
        alleles_df = format_data_for_optimization( # UPDATE THIS
            imputed_df)
        # Get loci that are not already a part of collection    
        loci_options = alleles_df.columns.difference(
            selected_loci.values_list('name', flat=True))
        alleles_df = alleles_df[list(loci_options)]
        # get form values
        randomize = form.cleaned_data['randomize']
        include_existing = form.cleaned_data['include_existing']
        metadata_cat = form.cleaned_data['metadata']
        n_targets = form.cleaned_data['n_targets']
        meta_df = project.get_metadata_table()[metadata_cat.name] if metadata_cat else None

        if include_existing:
            # Get starting resolution from exsiting targets
            starting_resolution = get_resolution(
                alleles_df[alleles_df.columns.intersection(
                    selected_loci.values_list('name', flat=True))])
        else:
            starting_resolution = get_resolution(get_starting_resolution(alleles_df.index))
        add_targets = optimization_loop(alleles_df, starting_resolution, n_targets, randomize, meta_df)
        add_target_inst = Locus.objects.filter(project=project, name__in=add_targets)
        target_collection.loci.add(*add_target_inst)
        return response
    




def profile_list_view(request, pk):
    target_col = TargetCollection.objects.get(id=pk)
    target_table = target_col.get_targets_table()
    project = target_col.project
    metadata_table = project.get_metadata_table()
    df = metadata_table.join(target_table)
    df.index.name = "GenomeID"
    html = df.reset_index().to_html(
        table_id="datatable",
        classes=['display', 'table', 'table-hover'],
        justify="left", index=False)
    
    metadata_idx = [i for i in range(1, metadata_table.shape[1] + 1)]
    allele_idx = [i for i in range(metadata_table.shape[1] + 1, df.shape[1] + 1)]
    

    return render(
                request, 'targets/profile_table.html',
                {
                    'table': html,
                    'targetcollection': target_col,
                    'meta_idx': metadata_idx,
                    'allele_idx': allele_idx
                }
                )
                

def get_resolution_figure_html(request, pk):
    """
    Get HTML for parcats figure.
    """
    print("HTML", request.POST, pk)
    target_col = TargetCollection.objects.get(id=pk)
    alleles = target_col.get_imputed_alleles()
    res = target_col.get_resolution_table()
    meta = target_col.project.get_metadata_table()

    mc = request.POST.get('metadata')
    print("MC", mc)
    if mc:
        metacat = MetadataCategory.objects.get(id=mc).name
    else:
        metacat = None
        
    html = draw_par_cats(
        res, meta, alleles, metacat,
        font_size=int(request.POST.get('fontsize', 12)),
        palette=request.POST.get('palette', 'plasma'))
    return HttpResponse(html)


def get_default_figure_html(pk):
    target_col = TargetCollection.objects.get(id=pk)
    alleles = target_col.get_imputed_alleles()
    res = target_col.get_resolution_table()
    meta = target_col.project.get_metadata_table()
        
    html = draw_par_cats(res, meta, alleles)
    return html


def resolution_view(request, pk):
    target_col = TargetCollection.objects.get(id=pk)
    form = ResolutionForm(instance=target_col)
    default_html = get_default_figure_html(pk)
    return render(
        request, 'targets/resolution_figure.html',
        {'form':form, 'pk': pk, 'targetcollection': target_col,
        'default_html': default_html})






# def add_target_collection(request, pk):
#     form = TargetCollectionForm(project=pk)
#     if request.method == 'POST':
#         form = TargetCollectionForm(request.POST)
#         if form.is_valid():
#             project = Project.objects.get(id=pk)
#             targets = TargetCollection(

#             )

#             return redirect(
#                 'lims:verify_new_samples',
#                 event_id=event, sample_type=sample_type)
#     return render(request, 'lims/samples_new.html', {'form':form})

# class GenomeListView(ListView):
#     template_name = "targets/genome_list.html"
#     model = Genome
#     def get_queryset(self, **kwargs):
#         new_context = Project.objects.get(pk=self.kwargs['pk']).get_active_genomes()
#         return new_context
    
#     def get_context_data(self, **kwargs):
#         context = super().get_context_data(**kwargs)
#         project = Project.objects.get(pk=self.kwargs['pk'])
#         context['project'] = project
#         return context


# # class LocusListView(ListView):
# #     template_name = "targets/locus_list.html"
# #     model = Locus
# #     def get_queryset(self, **kwargs):
# #         new_context = Project.objects.get(pk=self.kwargs['pk']).loci.all()
# #         return new_context
    
# #     def get_context_data(self, **kwargs):
# #         context = super().get_context_data(**kwargs)
# #         project = Project.objects.get(pk=self.kwargs['pk'])
# #         context['project'] = project
# #         return context
    


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
    

# ============== DATA ==================================
# def load_project_data(request, pk):
#     if request.method == 'POST':
#         form = DataUploadForm(request.POST, request.FILES)
#         if form.is_valid():
#             file_name = request.FILES['file']
#             project = Project.objects.get(pk=pk)
#             # Set all columns and row to inactive
#             project.deactivate_all_genomes()
#             project.deactivate_all_loci()
#             project.deactivate_all_metadata_categories()

#             organism = project.organism
            
#             df = pd.read_excel(file_name, index_col='id')
#             print(df.head())
#             # pull existing genome and loci data from database
#             # to filter out any values from the table that may
#             # not belong.
#             db_genomes = organism.get_genomes_list(request)
#             db_loci = organism.get_loci_list(request)
    
#             loci_cols = df.columns.intersection(db_loci)
#             genome_idx = df.index.intersection(db_genomes)

#             threshold_percentage = 75

#             # Calculate the minimum number of non-NaN values required based on the threshold percentage
#             threshold_count = int((100 - threshold_percentage) / 100 * len(loci_cols))

#             # Drop rows with more than the specified threshold count of NaN values
#             df = df.dropna(thresh=threshold_count)

#             # Consider all other columns to be metadata
#             metadata_cols = df.columns.difference(db_loci)
#             # keep only rows (genomes) that are found in the database
#             meta_df = df.loc[genome_idx, metadata_cols].to_dict()
#             alleles_df = df.loc[genome_idx, loci_cols].to_dict()
#             # TODO: add warning if alleles/genomes are removed
#             # add loci to database
#             loci_dict = {}
#             print("creating LOCI")
#             for locus in loci_cols:
#                 locus_inst, created = Locus.objects.get_or_create(
#                     name=locus,
#                     project=project,
#                 )
#                 if not created:
#                     # Created are active by default but if
#                     # locus already exists move to active
#                     locus_inst.active = True
#                     locus_inst.save()
#                 loci_dict[locus] = locus_inst
#             print("creating GENOMES")

#             # add genomes to database
#             genome_dict = {}
#             for genome in genome_idx:
#                 genome_inst, created = Genome.objects.get_or_create(
#                     name=genome,
#                     project=project,
#                 )
#                 if not created:
#                     # Created are active by default but if
#                     # genome already exists move to active
#                     genome_inst.active = True
#                     genome_inst.save()
#                 genome_dict[genome] = genome_inst
#             print("creating metacat")

#             # add metadata cateories to database
#             meta_dict = {}
#             for meta in metadata_cols:
#                 meta_inst, created = MetadataCategory.objects.get_or_create(
#                     name=meta,
#                     project=project,
#                 )
#                 if not created:
#                     # Created are active by default but if
#                     # metadata cat already exists move to active
#                     meta_inst.active = True
#                     meta_inst.save()
#                 meta_dict[meta] = meta_inst

#             # add alleles
#             print("creating alleles")

#             for locus, values in alleles_df.items():
#                 locus_inst = loci_dict[locus]
#                 for genome, allele in values.items():
#                     if str(allele) != 'nan': # don't add alleles with missing data
#                         Allele.objects.update_or_create(
#                             genome=genome_dict[genome],
#                             locus=locus_inst,
#                             project=project,
#                             defaults={'value': allele}
#                             )
#             print("creating METADATA")

#             # add metadata
#             for metacat, values in meta_df.items():
#                 metacat_inst = meta_dict[metacat]
#                 for genome, value in values.items():
#                     if str(value) != 'nan': # don't add metadata with missing values
#                         Metadata.objects.update_or_create(
#                             genome=genome_dict[genome],
#                             category=metacat_inst,
#                             project=project,
#                             defaults={'value': value}
#                             )
                       

#             return redirect('targets:project_detail', pk=pk)
#     else:
#         form = DataUploadForm()

#     return render(request, 'targets/data_upload.html', {'form': form})



# def data_list_view(request, pk):
#     project = Project.objects.get(id=pk)
#     active_genomes = project.get_active_genomes()
#     active_loci = project.get_active_loci()
#     active_metadata = project.get_active_metadata()
#     alleles = Allele.objects.filter(
#         project=project,
#         genome__in=active_genomes.values_list('id', flat=True),
#         locus__in=active_loci.values_list('id', flat=True)
#         ).values('genome', 'locus', 'value')
#     meta = Metadata.objects.filter(
#         project=project,
#         genome__in=active_genomes.values_list('id', flat=True),
#         category__in=active_metadata.values_list('id', flat=True)
#         ).values('genome', 'category', 'value')
#     allele_df = pd.DataFrame.from_dict(alleles)
#     allele_df = allele_df.pivot(index="genome", columns='locus', values='value')
#     idx_dict = {gen.id: gen.name for gen in active_genomes}
#     allele_df.index = [idx_dict.get(i) for i in allele_df.index]
#     col_dict = {loc.id: loc.name for loc in active_loci}
#     allele_df.columns = [col_dict.get(i) for i in allele_df.columns]

#     meta_df = pd.DataFrame.from_dict(meta)
#     meta_df = meta_df.pivot(index="genome", columns='category', values='value')
#     idx_dict = {gen.id: gen.name for gen in active_genomes}
#     meta_df.index = [idx_dict.get(i) for i in meta_df.index]
#     col_dict = {loc.id: loc.name for loc in active_metadata}
#     meta_df.columns = [col_dict.get(i) for i in meta_df.columns]

#     df = meta_df.join(allele_df)
#     df = df.fillna('')
#     print(df.head())
#     html = df.to_html(
#         table_id="datatable",
#         classes=['display', 'table', 'table-hover'],
#         justify="left", index=True)
#     print(html[:100])

#     return render(
#                 request, 'targets/data_table.html',
#                 {'table': html, 'project': project}
#                 )