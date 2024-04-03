
from .models import *
from django.forms import ModelForm, CheckboxInput, SelectMultiple
from django import forms
from django.utils.safestring import mark_safe
from django.core.validators import FileExtensionValidator
from django.forms.widgets import CheckboxSelectMultiple


class PrimerCollectionForm(ModelForm):
    class Meta:
        model = PrimerCollection
        fields = ['name', 'description', 'notes']
        
class AmpliconTargetForm(ModelForm):
    class Meta:
        model = AmpliconTarget
        fields = ['name', 'primer_collection', 'locus', ]