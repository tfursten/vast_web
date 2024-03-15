from django.db import models
from Bio.Seq import Seq
from targets.models import Locus, Project


class PrimerCollection(models.Model):
    """
    Set of amplicon targets and primers for a single reaction
    """
    name = models.CharField(max_length=255, blank=True, null=True, unique=True)
    description = models.CharField(max_length=255, blank=True, null=True)
    notes = models.TextField(null=True, blank=True)
    created_on = models.DateTimeField(auto_now_add=True, db_index=True, null=True)

    @property
    def count_loci(self):
        """
        Property returning the number of Amplicon Targets associated with the primer collection
        """
        return self.amplicontarget_set.filter().count()
    


class AmpliconTarget(models.Model):
    """
    Amplicon target region
    """
    name = models.CharField(max_length=255, blank=True, null=True)
    primer_collection = models.ForeignKey(PrimerCollection, on_delete=models.CASCADE)
    locus = models.ForeignKey(Locus, null=True, blank=False, on_delete=models.SET_NULL)
    target_start = models.PositiveIntegerField(blank=True, null=True)
    target_end = models.PositiveIntegerField(blank=True, null=True)
    forward_primer_zone_start = models.PositiveIntegerField(blank=True, null=True)
    forward_primer_zone_end = models.PositiveIntegerField(blank=True, null=True)
    reverse_primer_zone_start = models.PositiveIntegerField(blank=True, null=True)
    reverse_primer_zone_end = models.PositiveIntegerField(blank=True, null=True)



class Sequence(models.Model):
    """
    A sequence for an amplicon target allele
    """
    name = models.CharField(max_length=255, blank=True, null=True)
    seq = models.TextField()
    amp = models.ForeignKey(
        AmpliconTarget, on_delete=models.CASCADE)

class Primer(models.Model):
    """
    Primer object
    """
    name = models.CharField(max_length=255, blank=True, null=True)
    amp = models.ForeignKey(AmpliconTarget, blank=True, null=True, on_delete=models.CASCADE)
    seq = models.CharField(max_length=255, blank=True, null=True)
    flank = models.CharField(max_length=16, choices=[('Forward', 'Forward'), ('Reverse', 'Reverse')], blank=True, null=True)
    location = models.PositiveIntegerField(blank=True, null=True)
    user_added = models.BooleanField(default=False) # differentiate between manually added primers and optimized primers
    
    @property
    def length(self):
        """
        Get primer length
        """
        return len(self.seq)
    
    @property
    def reverse_complement(self):
        """
        Get reverse complement of primer sequence
        """
        seq = Seq(self.seq)
        return str(seq.reverse_complement())
    

    