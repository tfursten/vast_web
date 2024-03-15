# Generated by Django 5.0.2 on 2024-03-14 22:04

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('targets', '0002_locus_reference_locus_size_locus_start_position'),
    ]

    operations = [
        migrations.CreateModel(
            name='AmpliconTarget',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('target_start', models.PositiveIntegerField(blank=True, null=True)),
                ('target_end', models.PositiveIntegerField(blank=True, null=True)),
                ('forward_primer_zone_start', models.PositiveIntegerField(blank=True, null=True)),
                ('forward_primer_zone_end', models.PositiveIntegerField(blank=True, null=True)),
                ('reverse_primer_zone_start', models.PositiveIntegerField(blank=True, null=True)),
                ('reverse_primer_zone_end', models.PositiveIntegerField(blank=True, null=True)),
                ('locus', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='targets.locus')),
            ],
        ),
        migrations.CreateModel(
            name='Primer',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('seq', models.CharField(blank=True, max_length=255, null=True)),
                ('flank', models.CharField(blank=True, choices=[('Forward', 'Forward'), ('Reverse', 'Reverse')], max_length=16, null=True)),
                ('location', models.PositiveIntegerField(blank=True, null=True)),
                ('user_added', models.BooleanField(default=False)),
                ('amp', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='primers.amplicontarget')),
            ],
        ),
        migrations.CreateModel(
            name='PrimerCollection',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('description', models.CharField(blank=True, max_length=255, null=True)),
                ('notes', models.TextField(blank=True, null=True)),
                ('project', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='targets.project')),
            ],
        ),
        migrations.AddField(
            model_name='amplicontarget',
            name='primer_collection',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='primers.primercollection'),
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('seq', models.TextField()),
                ('amp', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='primers.amplicontarget')),
            ],
        ),
    ]