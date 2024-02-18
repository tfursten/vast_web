from io import StringIO
import requests
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
        if project.organism.id != organism.id:
            project.loci.clear()
            project.genomes.clear()
        # Save any updates to the database fields
        project.organism = organism
        project.save()
        if created:
            messages.success(
                request, "{0} sucessfully created and added to project {1}".format(
                    organism.name, project.name))
        else:
            messages.success(
                request, "{0} sucessfully updated and added to project {1}".format(
                    organism.name, project.name))
        return redirect('targets:project_detail', pk=pk)
        
    return render(request, "targets/select_organism.html")


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
    for genus in res:
        isolates = False
        seqdef = False
        species = False
        if genus['name'] == "rMLST":
            continue
        for db in genus['databases']:
    #         print(database['name'])
            if db['name'].split("_")[-1] == 'isolates':
                isolates = db['name']
                species =" ".join(db['description'].split(" ")[:-1])
            if db['name'].split("_")[-1] == "seqdef":
                seqdef = db['name']
        if isolates and seqdef:
            data.append({
                'id': idx, 
                'species': species,
                'isolates': isolates,
                'seqdef': seqdef})
            idx += 1
    return JsonResponse(data, safe=False)

# class OrganismUpdateView(UpdateView):
#     model = Panel
#     template_name_suffix = '_update'
#     form_class = PanelForm
#     success_message = "Panel was successfully updated:  %(name)s"
    
#     def get_success_url(self):
#         return reverse('targets:project_detail', args=(self.object.id,))



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
        'https://rest.pubmlst.org/db/{database}/loci?return_all=1'.format(
            database=project.organism.db_isolates), params=request.GET)
    split = str.split
    res = json.loads(r.text)
    data = [
        {'id': i, 'description': v} for i, v in enumerate(
        map(lambda x: split(x, "/")[-1], res['loci']))]

    return JsonResponse(data, safe=False)


def select_loci_list(request, pk):
    if request.method == 'POST':
        print(request.POST)
        loci_resp = request.POST
        created_count = 0
        updated_count = 0
        project = Project.objects.get(id=pk)
        project.loci.clear()
        for locus in loci_resp.getlist('loci'):
            print(locus)
            locus, created = Locus.objects.update_or_create(
                    name=locus,
                    organism=project.organism,
                )
            if created:
                created_count += 1
            else:
                updated_count += 1
            project.loci.add(locus)
        project.save()
        msg = "{0} loci sucessfully added to project {1}".format(created_count, project.name) if created_count else ""
        msg += "{0} successfully updated for project {1}".format(updated_count, project.name) if updated_count else ""
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
    
