from django.db import models

class Organism(models.Model):
    name = models.CharField(max_length=255)
    db_isolates = models.CharField(max_length=255)
    db_seqdef = models.CharField(max_length=255)

    @property
    def name_no_spaces(self):
        return self.name.replace("-like", "").replace(' ', '-').replace(".", "").replace("/", "")
    
    def __str__(self):
        return f"{self.name}"

class MetadataCategory(models.Model):
    name = models.CharField(max_length=255)
    def __str__(self):
        return f"{self.name}"


class Metadata(models.Model):
    category = models.ForeignKey(MetadataCategory, on_delete=models.CASCADE)
    value = models.CharField(max_length=255)

    def __str__(self):
        return f"{self.category} (Value: {self.value})"

class Locus(models.Model):
    name = models.CharField(max_length=255)
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE)

    @property
    def db_fasta_path(self):
        return 'https://rest.pubmlst.org/db/{0}/loci/{1}/alleles_fasta'.format(
            self.organism.db_seqdef, self.name)
    
    @property
    def count_alleles(self):
        return self.allele_set.count()
    
    def __str__(self):
        return f"{self.name} (Org: {self.organism.name})"

class Allele(models.Model):
    name = models.CharField(max_length=255)
    locus = models.ForeignKey(Locus, on_delete=models.CASCADE)
    
    def __str__(self):
        return f"Locus: {self.locus.name} (ID: {self.name})"

class Genome(models.Model):
    name = models.CharField(max_length=255)
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE)
    db_id = models.PositiveIntegerField()
    alleles = models.ManyToManyField(Allele)
    metadata = models.ManyToManyField(Metadata)

    def count_metadata(self):
        return self.metadata.count()
    
    def count_alleles(self):
        return self.alleles.count()
    
    def count_unique_loci(self):
        return Locus.objects.filter(allele__genome=self).distinct().count()
    
    def metadata_fields(self):
        return MetadataCategory.objects.filter(metadata__genome=self).distinct()
    
    def __str__(self):
        return f"{self.name} (ID: {self.db_id})"


class Project(models.Model):
    name = models.CharField(max_length=255)
    description = models.CharField(max_length=255, null=True, blank=True)
    notes = models.TextField(null=True, blank=True)
    loci = models.ManyToManyField(Locus)
    genomes = models.ManyToManyField(Genome)
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True, blank=True)
    created_on = models.DateTimeField(auto_now_add=True, db_index=True, null=True)


    @property
    def count_loci(self):
        return self.loci.count()
    
    @property
    def count_genomes(self):
        return self.genomes.count()
    
    @property
    def count_metadata(self):
        return MetadataCategory.objects.filter(metadata__genome__project=self).distinct().count()

    def __str__(self):
        return f"{self.name}"
    
    class Meta:
        ordering = ("-created_on",)


class TargetCollection(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    loci = models.ManyToManyField(Locus)
    created_on = models.DateTimeField(auto_now_add=True, db_index=True, null=True)


    @property
    def count_loci(self):
        return self.loci.count()
    
    def get_resolution(self):
        # TODO: calculate resolution of current loci given organism genotypes in project
        pass

    class Meta:
        ordering = ("-created_on",)









