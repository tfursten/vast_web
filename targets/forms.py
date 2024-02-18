from django import forms
from django.forms import ModelForm, CheckboxInput, SelectMultiple
from .models import *



class ProjectForm(ModelForm):
    class Meta:
        model = Project
        fields = ['name', 'description', 'notes']
        widgets = {

        }
    
    
  