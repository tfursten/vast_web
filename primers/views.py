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


