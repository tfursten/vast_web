# Generated by Django 5.0.2 on 2024-02-22 04:26

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('targets', '0006_project_data_delete_inputdata'),
    ]

    operations = [
        migrations.AddField(
            model_name='locus',
            name='alleles',
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
    ]
