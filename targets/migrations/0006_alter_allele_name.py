# Generated by Django 5.0.2 on 2024-02-20 07:15

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('targets', '0005_alter_metadatacategory_name'),
    ]

    operations = [
        migrations.AlterField(
            model_name='allele',
            name='name',
            field=models.PositiveIntegerField(),
        ),
    ]
