from django import forms
import pandas as pd
from django.forms import ModelForm, CheckboxInput, SelectMultiple
from .models import *



class ProjectForm(ModelForm):
    class Meta:
        model = Project
        fields = ['name', 'description', 'notes']
        widgets = {

        }
    
    
class DataUploadForm(forms.Form):
    file = forms.FileField(label='Upload BIGSdb excel file from PubMLST')

    def clean(self):
        uploaded_file = self.cleaned_data.get('file')
        if uploaded_file:
            try:
                pd.read_excel(uploaded_file, nrows=1, usecols=[0])
            # Validate that the uploaded file is a CSV file
            except Exception as exc:
                raise forms.ValidationError('File must be an Excel file.') from exc
        return uploaded_file
    


