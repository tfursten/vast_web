from django import forms
import pandas as pd
from django.forms import ModelForm, CheckboxInput, SelectMultiple
from .models import *
from django.core.validators import FileExtensionValidator
from django.forms.widgets import CheckboxSelectMultiple
from django.utils.safestring import mark_safe
import plotly.express as px
from django.db.models import Q



class ProjectForm(ModelForm):
    class Meta:
        model = Project
        fields = ['name', 'description', 'notes']
        widgets = {

        }

class TargetCollectionForm(ModelForm):
    class Meta:
        model = TargetCollection
        fields = ['name', 'project', 'description']
    
        labels = {
            'project': mark_safe("""
                        <b>Select Project </b>
                        <a href="#" data-toggle="tooltip"
                        title='Only projects that have an assigned organism and data uploaded are displayed. If project is not in list, finish project setup before adding target collection.'>
                        <span data-feather="info"
                        style="display: inline-block;"></span></a>"""),
            
        }
        
    def __init__(self, *args, **kwargs):
        if 'project' in kwargs:
            project = kwargs.pop('project')
            kwargs.update(initial={
                'project': project
            })
        super(TargetCollectionForm, self).__init__(*args, **kwargs)
        self.fields['project'].queryset = Project.objects.filter(
            id__in=Locus.objects.filter(active=True).values_list('project', flat=True).distinct())
        

  
    
            



class TargetCollectionUpdateForm(ModelForm):
    """
    Update form removes project. Once a target collection
    is created, don't allow project to be changed.
    """
    class Meta:
        model = TargetCollection
        fields = ['name', 'description']

    
class DataUploadForm(forms.Form):
    file = forms.FileField(label='Upload BIGSdb file from PubMLST')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Add a FileExtensionValidator to the file_field
        self.fields['file'].validators.append(FileExtensionValidator(allowed_extensions=['xlsx', 'csv']))


class AddTargets(ModelForm):
    class Meta:
        model = TargetCollection
        exclude = ['name', 'project', 'description', 'loci']

    n_targets = forms.IntegerField(min_value=1, label=mark_safe("<b>Number of targets to add</b>"), initial=1)
    metadata = forms.ModelChoiceField(
        label=mark_safe("""
                        <b>Differentiate genomes by metadata category </b>
                        <a href="#" data-toggle="tooltip"
                        title='Targets are selected that best differentiate genomes by selected feature. Genomes with missing metadata are not included in the calculation.'>
                        <span data-feather="info"
                        style="display: inline-block;"></span></a>"""),
        queryset=MetadataCategory.objects.filter(active=True), widget=forms.RadioSelect, required=False)
    include_existing = forms.BooleanField(
        label=mark_safe("""
                        <b>Include existing targets in calculation </b>
                        <a href="#" data-toggle="tooltip"
                        title='New targets will be calculated based on resolution with any existing targets. Ignoring existing targets will start with no resolution but will not add any targets that are currently in the collection. The later option is good for building in redundancy.'>
                        <span data-feather="info"
                        style="display: inline-block;"></span></a>"""),
                initial=True, required=False)
    randomize = forms.BooleanField(
        label=mark_safe("""
                        <b>Randomize results </b>
                        <a href="#" data-toggle="tooltip"
                        title='In the case of a tie, a target will be picked at random instead of deterministically.'>
                        <span data-feather="info"
                        style="display: inline-block;"></span></a>"""),

                initial=False, required=False)


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        instance = kwargs.get('instance')

        if instance:
            # Filter categories based on the primary key
            print('inst', instance)
            self.fields['metadata'].queryset = MetadataCategory.objects.filter(
                active=True, selected=True, project=instance.project)
        else:
            # If no instance is provided, show all categories
            self.fields['metadata'].queryset = MetadataCategory.objects.filter(
                active=True)
    


class ResolutionForm(ModelForm):
    class Meta:
        model = TargetCollection
        exclude = ['name', 'project', 'description', 'loci']

    metadata = forms.ModelChoiceField(
        label="Color by metadata category",
        queryset=MetadataCategory.objects.filter(active=True), required=False)
    fontsize = forms.IntegerField(min_value=8, max_value=20, initial=12)
    palette = forms.ChoiceField(choices=[(n, n.capitalize()) for n in px.colors.named_colorscales()], initial='plasma')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        instance = kwargs.get('instance')

        if instance:
            # Filter categories based on the primary key
            self.fields['metadata'].queryset = MetadataCategory.objects.filter(
                active=True, project=instance.project)
        else:
            # If no instance is provided, show all categories
            self.fields['metadata'].queryset = MetadataCategory.objects.filter(
                active=True)
    
