from django.db import models
import requests
import json
import pandas as pd
from .utils import get_resolution, format_data_for_optimization

class Organism(models.Model):
    """
    Model representing an PubMLST organism and path to databases in API
    """
    name = models.CharField(max_length=255)
    db_isolates = models.CharField(max_length=255)
    db_seqdef = models.CharField(max_length=255)

    @property
    def name_no_spaces(self):
        """
        Return organism name without spaces for adding to paths
        """
        return self.name.replace("-like", "").replace(' ', '-').replace(".", "").replace("/", "")
    
    def __str__(self):
        """
        String representation of the organism.
        """
        return f"{self.name}"
    
    def get_genomes_list(self, request):
        """
        Pull data on all genomes for organism from PubMLST API
        """
        r = requests.get(
        f'https://rest.pubmlst.org/db/{self.db_isolates}/genomes?return_all=1', params=request.GET)
        split = str.split
        res = json.loads(r.text)
        return list(map(lambda x: int(split(x, "/")[-1]), res['isolates']))

    def get_loci_list(self, request):
        """
        Pull data on all loci for organism from PubMLST API.
        """
        r = requests.get(
            f'https://rest.pubmlst.org/db/{self.db_isolates}/loci?return_all=1', params=request.GET)
        split = str.split
        res = json.loads(r.text)
        return list(map(lambda x: split(x, "/")[-1], res['loci']))




class Project(models.Model):
    """
    Model representing a project.
    """
    name = models.CharField(max_length=64)
    description = models.CharField(max_length=255, null=True, blank=True)
    notes = models.TextField(null=True, blank=True)
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True, blank=True)
    data = models.FileField(upload_to='uploads/', null=True, blank=True)
    imputed_data = models.FileField(upload_to='uploads/', null=True, blank=True)
    created_on = models.DateTimeField(auto_now_add=True, db_index=True, null=True)
    
    class Meta:
        """
        Metadata for ordering projects by creation date.
        """
        ordering = ("-created_on",)

    def __str__(self):
        """
        String representation of the project.
        """
        return f"{self.name}"

    @property
    def count_loci(self):
        """
        Property returning the count of loci associated with the project.
        """
        return self.locus_set.count()

    @property
    def count_genomes(self):
        """
        Property returning the count of genomes associated with the project.
        """
        return self.genome_set.count()
    
    @property
    def count_metadata(self):
        """
        Property returning the count of unique metadata categories associated with the project.
        """
        return self.metadatacategory_set.count()

    @property
    def count_active_loci(self):
        """
        Property returning the count of loci associated with the project.
        """
        return self.locus_set.filter(active=True).count()

    @property
    def count_active_genomes(self):
        """
        Property returning the count of genomes associated with the project.
        """
        return self.genome_set.filter(active=True).count()
    
    @property
    def count_active_metadata(self):
        """
        Property returning the count of unique metadata categories associated with the project.
        """
        return self.metadatacategory_set.filter(active=True).count()
    
    @property
    def count_target_collections(self):
        """
        Property returning the number of target collections associated with the project.
        """
        return self.targetcollection_set.filter().count()
    
    def get_genomes(self):
        """
        Return genomes associated with the project
        """
        return self.genome_set.all()

    def get_loci(self):
        """
        Return loci associated with the project
        """
        return self.locus_set.all()

    def get_metadata(self):
        """
        Return metadata categories associated with the project.
        """
        return self.metadatacategory_set.all()
    
    def get_active_genomes(self):
        """
        Return genomes associated with the project
        """
        return self.genome_set.filter(active=True)

    def get_active_loci(self):
        """
        Return loci associated with the project
        """
        return self.locus_set.filter(active=True)

    def get_active_metadata(self):
        """
        Return metadata categories associated with the project.
        """
        return self.metadatacategory_set.filter(active=True)
    
    def get_deactive_genomes(self):
        """
        Return genomes associated with the project
        """
        return self.genome_set.filter(active=False)

    def get_deactive_loci(self):
        """
        Return loci associated with the project
        """
        return self.locus_set.filter(active=False)

    def get_deactive_metadata(self):
        """
        Return metadata categories associated with the project.
        """
        return self.metadatacategory_set.filter(active=False)   
     
    def deactivate_all_genomes(self):
        """
        Deactivate all Genome instances associated with this Project.
        """
        Genome.objects.filter(project=self, active=True).update(active=False)


    def deactivate_all_loci(self):
        """
        Deactivate all Locus instances associated with this Project.
        """
        Locus.objects.filter(project=self, active=True).update(active=False)
   

    def deactivate_all_metadata_categories(self):
        """
        Deactivate all metadata category instances associated with this Project.
        """
        MetadataCategory.objects.filter(project=self, active=True).update(active=False)

    def get_allele_table(self):
        """
        return dataframe using active genomes and loci
        """
        active_genomes = list(self.get_active_genomes().values_list('name', flat=True))
        active_loci = list(self.get_active_loci().values_list('name', flat=True))
        df = pd.read_excel(self.data, index_col=0, na_values='', keep_default_na=False, dtype=str)
        return df.loc[active_genomes, active_loci].fillna('')

    def get_metadata_table(self):
        """
        return dataframe using active genomes and metadata
        """
        active_genomes = list(self.get_active_genomes().values_list('name', flat=True))
        active_metadata = list(self.get_active_metadata().values_list('name', flat=True))
        df = pd.read_excel(self.data, index_col=0, na_values='', keep_default_na=False, dtype=str)
        return df.loc[active_genomes, active_metadata].fillna('')
    
    def get_imputed_data(self):
        """
        Return dateframe of imputed values. The multindex is a combination
        of the genome id and the locus id and the value is the imputed allele.
        """
        return pd.read_csv(self.imputed_data, index_col=[0, 1], sep="\t")

class MetadataCategory(models.Model):
    """
    Model representing a metadata category.
    Represents a column in metadata table.
    Active state determines whether column is added to
    table.
    """
    name = models.CharField(max_length=255)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    active = models.BooleanField(default=True)
    genomes = models.PositiveIntegerField(null=True, blank=True) # number of genomes with metadata
    values = models.PositiveIntegerField(null=True, blank=True) # number of unique values

    class Meta:
        unique_together = ('project', 'name')

    def __str__(self):
        """
        String representation of the metadata category.
        """
        return f"{self.name}"



class Locus(models.Model):
    """
    Model representing a locus.
    Represents a column in allele table.
    Active state determines whether column is added to
    table.
    """
    name = models.CharField(max_length=255)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    active = models.BooleanField(default=True)
    alleles = models.PositiveIntegerField(null=True, blank=True) # number of unqiue alleles

    class Meta:
        unique_together = ('project', 'name')

    @property
    def db_fasta_path(self):
        """
        Property returning the URL for the locus's alleles in FASTA format.
        """
        return 'https://rest.pubmlst.org/db/{0}/loci/{1}/alleles_fasta'.format(
            self.organism.db_seqdef, self.name)
    
    # @property
    # def count_alleles(self):
    #     """
    #     Property returning the count of unique alleles associated with the locus.
    #     """
    #     return self.allele_set.all().values_list('value', flat=True).distinct().count()
    
    def __str__(self):
        """
        String representation of the locus.
        """
        return f"{self.name} (Org: {self.project.organism.name})"

class Genome(models.Model):
    """
    Model representing a genome.
    Represents index in allele and metadata tables.
    Active state determines whether row is added to
    table.
    """
    name = models.PositiveIntegerField()
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    active = models.BooleanField(default=True)
    alleles = models.PositiveIntegerField(null=True, blank=True) # number of loci with alleles
    metadata = models.PositiveIntegerField(null=True, blank=True) # number of metadata categories with data
    

    # def count_metadata(self):
    #     """
    #     Method returning the count of metadata associated with the genome.
    #     """
    #     return self.metadata_set.count()
    
    # def count_alleles(self):
    #     """
    #     Method returning the count of alleles associated with the genome.
    #     """
    #     return self.allele_set.count()
    
    def __str__(self):
        """
        String representation of the genome.
        """
        return f"{self.name}"

    class Meta:
        ordering = ("name",)
        unique_together = ('project', 'name')

class TargetCollection(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    description = models.TextField(max_length=255, blank=True, null=True)
    loci = models.ManyToManyField(Locus)
    created_on = models.DateTimeField(auto_now_add=True, db_index=True, null=True)


    @property
    def count_loci(self):
        return self.loci.count()
   
    def get_resolution_table(self):
        project = self.project
        alleles_df = project.get_allele_table()
        imputed_df = project.get_imputed_data()
        selected_loci = self.loci.values_list('name', flat=True)
        alleles_df = format_data_for_optimization(
            alleles_df, imputed_df)[selected_loci]
        res = []
        for i in range(alleles_df.shape[1]):
            print(i)
            print(alleles_df.iloc[:, range(0, i+1)])
            res.append(get_resolution(alleles_df.iloc[:, range(0, i+1)]))
        return pd.DataFrame(res, index=alleles_df.columns, columns=alleles_df.index)


    def get_targets_table(self):
        project = self.project
        active_genomes = list(project.get_active_genomes().values_list('name', flat=True))
        loci = list(self.loci.all().values_list('name', flat=True))
        df = pd.read_excel(project.data, index_col=0, na_values='', keep_default_na=False, dtype=str)
        return df.loc[active_genomes, loci].fillna('')


    class Meta:
        ordering = ("-created_on",)
        unique_together = ('project', 'name')






# class Allele(models.Model):
#     """
#     Model representing an allele.
#     Represents values of allele table. 
#     """
#     locus = models.ForeignKey(Locus, on_delete=models.CASCADE)
#     genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
#     project = models.ForeignKey(Project, on_delete=models.CASCADE)
#     value = models.CharField(max_length=64)

#     class Meta:
#         unique_together = ('project', 'locus', 'genome')
#         indexes = [
#             models.Index(fields=['project', 'locus', 'genome'])
#         ]
    
#     def __str__(self):
#         """
#         String representation of the allele.
#         """
#         return f"Locus: {self.locus.name} (ID: {self.value} Genome {self.genome.name} Project {self.project.id})"


# class Metadata(models.Model):
#     """
#     Model representing metadata.
#     Represents values in metadata table.
#     """
#     category = models.ForeignKey(MetadataCategory, on_delete=models.CASCADE)
#     genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
#     project = models.ForeignKey(Project, on_delete=models.CASCADE)
#     value = models.CharField(max_length=255)

#     class Meta:
#         unique_together = ('project', 'category', 'genome')
#         indexes = [
#             models.Index(fields=['project', 'category', 'genome'])
#         ]

#     def __str__(self):
#         """
#         String representation of the metadata.
#         """
#         return f"{self.category} (Value: {self.value}, Genome {self.genome.name} Project {self.project.id})"











