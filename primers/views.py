from django.shortcuts import render
from .models import *
from .forms import *

from django.views.generic import (
    ListView, CreateView, DeleteView, UpdateView, DetailView)
from django.contrib.messages.views import SuccessMessageMixin
from django.urls import reverse, reverse_lazy

# ============== PRIMER SET ====================

class PrimerCollectionListView(ListView):
    template_name_suffix = "_list"
    context_object_name = 'primercollection_list'
    model = PrimerCollection


class PrimerCollectionFormView(SuccessMessageMixin, CreateView):
    model = PrimerCollection
    template_name_suffix = '_new'
    form_class = PrimerCollectionForm
    success_message = "Primer collection was successfully added: %(name)s"

    def get_success_url(self):
        return reverse('primers:primer_collection_detail', args=(self.object.id,))


class PrimerCollectionDetailView(DetailView):
    model = PrimerCollection

class PrimerCollectionUpdateView(SuccessMessageMixin, UpdateView):
    model = PrimerCollection
    template_name_suffix = '_update'
    form_class = PrimerCollectionForm
    success_message = "Primer collection was successfully updated:  %(name)s"
    
    def get_success_url(self):
        return reverse('primers:primer_collection_detail', args=(self.object.id,))
    
    def get_initial(self):
        initial = super().get_initial()
        # Retrieve the object being updated
        obj = self.get_object()
        # Populate the initial data with values from the object
        initial['name'] = obj.name
        initial['description'] = obj.description
        initial['notes'] = obj.notes
        return initial


class PrimerCollectionDeleteView(DeleteView):
    model = PrimerCollection
    success_url = reverse_lazy('primers:primer_collection_list', args=())
    
    def post(self, request, *args, **kwargs):
        try:
            return self.delete(request, *args, **kwargs)
        except ProtectedError:
            return render(request, "primers/protected_error.html")


# ================== AMPLICON TARGETS ===============================

class AmpliconTargetListView(ListView):
    template_name_suffix = "_list"
    context_object_name = 'amplicontarget_list'
    model = AmpliconTarget


class AmpliconTargetFormView(SuccessMessageMixin, CreateView):
    model = AmpliconTarget
    template_name_suffix = '_new'
    form_class = AmpliconTargetForm
    success_message = "Amplicon target was successfully added: %(name)s"

    def get_success_url(self):
        return reverse('primers:amplicon_target_detail', args=(self.object.id,))


def sample_storage_label_options(request):
    form = SamplePrint(initial={
        'label_paper': 3, 'abbreviate': True,
        'sort_by1': 'LOCATION', 'sort_by2': 'EVENT', 'sort_by3': 'GRADE', 'sort_by4': 'NAME' })
    samples = Sample.objects.filter(collection_status="Collected").values(
        'id', 'name', 'subject__subject_ui', 'subject__first_name',
        'subject__last_name', 'subject__grade',
        'collection_event__name', 'location__name',
        'sample_type'
        )
    context = {'data': json.dumps(list(samples)), 'form': form}
    if request.method == "POST":
        form = SamplePrint(request.POST)
        if form.is_valid():
            start_position = request.POST.get('start_position')
            label_paper = request.POST.get('label_paper')
            abbreviate = request.POST.get('abbreviate')
            reps = request.POST.get('replicates')
            sort_by1 = request.POST.get('sort_by1')
            sort_by2 = request.POST.get('sort_by2')
            sort_by3 = request.POST.get('sort_by3')
            sort_by4 = request.POST.get('sort_by4')
            selected_samples = request.POST.getlist('ids')
            return sample_labels_pdf(samples=selected_samples,
                start_position=start_position,
                label_paper=label_paper, abbreviate=abbreviate, replicates=reps,
                sort_by1=sort_by1, sort_by2=sort_by2, sort_by3=sort_by3,
                sort_by4=sort_by4, outfile_name="storage_labels")
    return render(request, 'lims/sample_storage_print_options.html', context=context)   



class AmpliconTargetDetailView(DetailView):
    model = AmpliconTarget

class AmpliconTargetUpdateView(SuccessMessageMixin, UpdateView):
    model = AmpliconTarget
    template_name_suffix = '_update'
    form_class = AmpliconTargetForm
    success_message = "Amplicon target was successfully updated:  %(name)s"
    
    def get_success_url(self):
        return reverse('primers:amplicon_target_detail', args=(self.object.id,))
    
    def get_initial(self):
        initial = super().get_initial()
        # Retrieve the object being updated
        obj = self.get_object()
        # Populate the initial data with values from the object
        initial['name'] = obj.name
        initial['description'] = obj.description
        initial['notes'] = obj.notes
        return initial


class AmpliconTargetDeleteView(DeleteView):
    model = AmpliconTarget
    success_url = reverse_lazy('primers:primer_collection_list', args=())
    
    def post(self, request, *args, **kwargs):
        try:
            return self.delete(request, *args, **kwargs)
        except ProtectedError:
            return render(request, "primers/protected_error.html")


