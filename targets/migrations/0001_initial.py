# Generated by Django 3.2.25 on 2024-03-06 22:43

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=64)),
                ('description', models.CharField(blank=True, max_length=255, null=True)),
                ('notes', models.TextField(blank=True, null=True)),
                ('organism', models.CharField(blank=True, max_length=64, null=True)),
                ('created_on', models.DateTimeField(auto_now_add=True, db_index=True, null=True)),
            ],
            options={
                'ordering': ('-created_on',),
            },
        ),
        migrations.CreateModel(
            name='MetadataCategory',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('active', models.BooleanField(default=True)),
                ('selected', models.BooleanField(default=True)),
                ('n_genomes', models.PositiveIntegerField(blank=True, null=True)),
                ('n_values', models.PositiveIntegerField(blank=True, null=True)),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'ordering': ('name',),
            },
        ),
        migrations.CreateModel(
            name='Locus',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('active', models.BooleanField(default=True)),
                ('selected', models.BooleanField(default=True)),
                ('n_alleles', models.PositiveIntegerField(blank=True, null=True)),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'unique_together': {('project', 'name')},
            },
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.PositiveIntegerField()),
                ('active', models.BooleanField(default=True)),
                ('selected', models.BooleanField(default=True)),
                ('n_alleles', models.PositiveIntegerField(blank=True, null=True)),
                ('n_metadata', models.PositiveIntegerField(blank=True, null=True)),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'ordering': ('name',),
                'unique_together': {('project', 'name')},
            },
        ),
        migrations.CreateModel(
            name='TargetCollection',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('description', models.TextField(blank=True, max_length=255, null=True)),
                ('created_on', models.DateTimeField(auto_now_add=True, db_index=True, null=True)),
                ('loci', models.ManyToManyField(to='targets.Locus')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'ordering': ('-created_on',),
                'unique_together': {('project', 'name')},
            },
        ),
        migrations.CreateModel(
            name='Metadata',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.CharField(blank=True, max_length=16, null=True)),
                ('category', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.metadatacategory')),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.genome')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'unique_together': {('project', 'genome', 'category')},
            },
        ),
        migrations.CreateModel(
            name='Allele',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('allele', models.CharField(blank=True, max_length=16, null=True)),
                ('imputed', models.CharField(blank=True, max_length=16, null=True)),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.genome')),
                ('locus', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.locus')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='targets.project')),
            ],
            options={
                'unique_together': {('project', 'genome', 'locus')},
            },
        ),
    ]
