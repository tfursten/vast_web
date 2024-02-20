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

from django.contrib import messages
from django.db.models import ProtectedError, Count, Q, Max, Min
from .models import *
from .forms import *
from django.views.generic import (
    ListView, CreateView, DeleteView, UpdateView, DetailView)
from Bio import SeqIO, Align


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
def get_genomes_list(request, isolates_db):
    r = requests.get(
        f'https://rest.pubmlst.org/db/{isolates_db}/genomes?return_all=1', params=request.GET)
    split = str.split
    res = json.loads(r.text)
    return list(map(lambda x: int(split(x, "/")[-1]), res['isolates']))


def get_loci_list(request, isolates_db):
    r = requests.get(
        f'https://rest.pubmlst.org/db/{isolates_db}/loci?return_all=1', params=request.GET)
    split = str.split
    res = json.loads(r.text)
    return list(map(lambda x: split(x, "/")[-1], res['loci']))


def load_project_data(request, pk):
    if request.method == 'POST':
        form = DataUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file_name = request.FILES['file']
            project = Project.objects.get(pk=pk)
            project.loci.clear()
            project.genomes.clear()
            org = project.organism
            df = pd.read_excel(file_name, index_col='id')
            db_genomes = get_genomes_list(request, org.db_isolates)
            db_loci = get_loci_list(request, org.db_isolates)
            loci_cols = df.columns.intersection(db_loci)
            metadata_cols = df.columns.difference(db_loci)
            genome_idx = df.index.intersection(db_genomes)
            # keep only rows (genomes) that are found in the database
            df = df.loc[genome_idx]
            # add loci to database
            loci_allele_dict = {}
            add_loci = []
            for locus in loci_cols:
                locus_inst, _ = Locus.objects.get_or_create(
                    name=locus,
                    organism=project.organism,
                )
                add_loci.append(locus_inst)
                loci_allele_dict[locus] = {}
                if df[locus].dtype == object: # replace multiple alleles with only first one
                    df[locus] = df[locus].apply(lambda x: x.split(";")[0] if isinstance(x, str) else x)
                for allele in df[locus].dropna().unique():
                    if allele: # not none or nan
                        allele_inst, _ = Allele.objects.get_or_create(
                            name=allele,
                            locus=locus_inst
                        )
                        loci_allele_dict[locus][allele] = allele_inst

            # add metadata to database
            metadata_dict = {}
            for meta in metadata_cols:
                meta_cat_inst, _ = MetadataCategory.objects.get_or_create(
                    name=meta
                )

                metadata_dict[meta] = {}
             
                for value in df[meta].dropna().unique():
                    if value: # not none or nan
                        meta_inst, _ = Metadata.objects.get_or_create(
                            value=value,
                            category=meta_cat_inst
                        )
                        metadata_dict[meta][value] = meta_inst
            meta_df = df[metadata_cols]
            alleles_df = df[loci_cols]
            # create genomes
            add_genomes = []
            for genome in genome_idx:
                if genome: # not none or nan
                    genome_inst, _ = Genome.objects.get_or_create(
                        name=genome,
                        organism=project.organism,
                        )

                    # add metadata Many2Many relationships
                    genome_inst.alleles.clear() # start from scratch on new data upload
                    genome_inst.metadata.clear()
                    add_metadata = []
                    this_meta = meta_df.loc[genome].dropna()
                    for k, v in this_meta.items():
                        if v:
                            add_metadata.append(metadata_dict[k][v])
                    genome_inst.metadata.add(*add_metadata)
                    # add alleles Many2Many relationships
                    add_alleles = []
                    this_alleles = alleles_df.loc[genome].dropna()
                    for k, v in this_alleles.items():
                        if v:
                            add_alleles.append(loci_allele_dict[k][v])
                    genome_inst.alleles.add(*add_alleles)
                    genome_inst.save()
                    add_genomes.append(genome_inst)
            project.genomes.add(*add_genomes)
            project.loci.add(*add_loci)
            project.save()
            return redirect('targets:project_list')
    else:
        form = DataUploadForm()

    return render(request, 'targets/data_upload.html', {'form': form})



# ============== GENOME ==================================
class GenomeListView(ListView):
    template_name = "targets/genome_list.html"
    model = Genome
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).genomes.all()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    

def get_genomes_json(request, pk):
    project = Project.objects.get(id=pk)
    r = requests.get(
        f'https://rest.pubmlst.org/db/{project.organism.db_isolates}/genomes?return_all=1', params=request.GET)
    split = str.split
    res = json.loads(r.text)
    data = [
        {'id': i, 'description': v} for i, v in enumerate(
        map(lambda x: split(x, "/")[-1], res['isolates']))]

    return JsonResponse(data, safe=False)


def select_genomes_list(request, pk):
    if request.method == 'POST':
        genome_resp = request.POST
        created_count = 0
        project = Project.objects.get(id=pk)
        project.genomes.clear()
        new_genomes = []
        for g in genome_resp.getlist('genome'):
            if not Genome.objects.filter(name = g, organism=project.organism).exists():
                new_genomes.append(g)
            else:
                genome = Genome.objects.get(name=g, organism=project.organism)
                project.genomes.add(genome)
        old_genomes = len(genome_resp.getlist('genome')) - len(new_genomes)
        allele_urls = {k: f'https://rest.pubmlst.org/db/{project.organism.db_isolates}/isolates/{k}/allele_ids?return_all=1' for k in new_genomes }

        max_threads = 4  # PUBMLST REST API MAX
        def api_call(url):
            response = requests.get(url, params=request.GET)
            result = json.loads(response.text)['allele_ids']
            return result
        # Use ThreadPoolExecutor to make concurrent API calls
        with concurrent.futures.ThreadPoolExecutor(max_threads) as executor:
            # Use the map function to apply the api_call function to each URL in parallel
            results = dict(zip(allele_urls.keys(), list(executor.map(api_call, allele_urls.values()))))
        
        # Get dictionary of all unique loci and alleles to create/update
        uni_loci_alleles = {}
        for genome, result in results.items():
            for d in result:
                for locus, allele in d.items():
                    if isinstance(allele, list):
                        # sometimes the alleles are a list so pick first one
                        allele = allele[0]
                    if locus in uni_loci_alleles:
                        uni_loci_alleles[locus].add(allele)
                    else:
                        uni_loci_alleles[locus] = {allele}
        allele_objs = {}       
        # update/create loci and alleles
        for locus, alleles in uni_loci_alleles.items():
            loc, _ = Locus.objects.get_or_create(
                        name=locus,
                        organism=project.organism
                    )
            for allele in alleles:
                al, _ = Allele.objects.update_or_create(
                        name=allele,
                        locus=loc
                    )
                if locus in allele_objs:
                    allele_objs[locus][allele] = al
                else:
                    allele_objs[locus] = {allele: al}
        # Assign alleles to genomes
        add_genomes = []
        for genome, result in results.items():
            print(f"ON GENOME {genome}")
            gen, created = Genome.objects.get_or_create(
                    name=genome,
                    organism=project.organism,
                )
            add_alleles = []
            for d in result:
                for locus, allele in d.items():
                    if isinstance(allele, list):
                        # sometimes the alleles are a list so pick first one
                        allele = allele[0]
                    add_alleles.append(allele_objs[locus][allele])
            gen.alleles.add(*add_alleles)
            gen.save()
            if created:
                add_genomes.append(gen)
                created_count += 1
        # add updated genomes
        project.genomes.add(*add_genomes)
        project.save()

        msg = f"Successfully created {created_count} genomes and added to project {project.name}."if created_count else ""
        msg += f" Successfully updated {old_genomes} genomes for project {project.name}." if old_genomes else ""
        messages.success(request, msg)
        return redirect('targets:project_detail', pk=pk) # TODO: change redirect to genome list when created
        
    return render(request, "targets/select_genomes.html", {'pk': pk})

# ============== LOCI ==================================

class LocusListView(ListView):
    template_name = "targets/locus_list.html"
    model = Locus
    def get_queryset(self, **kwargs):
        new_context = Project.objects.get(pk=self.kwargs['pk']).loci.all()
        return new_context
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs['pk'])
        context['project'] = project
        return context
    


def get_loci_json(request, pk):
    project = Project.objects.get(id=pk)
    r = requests.get(
        f'https://rest.pubmlst.org/db/{project.organism.db_isolates}/loci?return_all=1', params=request.GET)
    split = str.split
    res = json.loads(r.text)
    data = [
        {'id': i, 'description': v} for i, v in enumerate(
        map(lambda x: split(x, "/")[-1], res['loci']))]

    return JsonResponse(data, safe=False)


def select_loci_list(request, pk):
    if request.method == 'POST':
        loci_resp = request.POST
        created_count = 0
        updated_count = 0
        project = Project.objects.get(id=pk)
        project.loci.clear()
        add_loci = []
        for locus in loci_resp.getlist('loci'):
            locus, created = Locus.objects.get_or_create(
                    name=locus,
                    organism=project.organism,
                )
            add_loci.append(locus)
            if created:
                created_count += 1
            else:
                updated_count += 1
        project.loci.add(*add_loci)
        project.save()
        msg = f"Successfully added {created_count} loci project {project.name}. " if created_count else ""
        msg += f"Successfully updated {updated_count} loci for project {project.name}" if updated_count else ""
        messages.success(request, msg)
        return redirect('targets:project_detail', pk=pk) # TODO: change redirect to loci list when created
        
    return render(request, "targets/select_loci.html", {'pk': pk})
    # r = requests.get('https://rest.pubmlst.org/db', params=request.GET)


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
    
