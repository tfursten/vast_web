from django import forms
import pandas as pd
from django.forms import ModelForm, CheckboxInput, SelectMultiple
from .models import *
from django.core.validators import FileExtensionValidator



class ProjectForm(ModelForm):
    class Meta:
        model = Project
        fields = ['name', 'description', 'notes']
        widgets = {

        }
    
    
class DataUploadForm(forms.Form):
    file = forms.FileField(label='Upload BIGSdb excel file from PubMLST')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Add a FileExtensionValidator to the file_field
        self.fields['file'].validators.append(FileExtensionValidator(allowed_extensions=['xlsx']))
