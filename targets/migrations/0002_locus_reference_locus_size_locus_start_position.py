# Generated by Django 5.0.2 on 2024-03-07 18:08

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('targets', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='locus',
            name='reference',
            field=models.CharField(blank=True, max_length=64, null=True),
        ),
        migrations.AddField(
            model_name='locus',
            name='size',
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='locus',
            name='start_position',
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
    ]
