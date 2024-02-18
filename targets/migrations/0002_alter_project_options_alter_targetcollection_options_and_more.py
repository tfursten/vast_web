# Generated by Django 5.0.2 on 2024-02-17 16:39

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('targets', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='project',
            options={'ordering': ('-created_on',)},
        ),
        migrations.AlterModelOptions(
            name='targetcollection',
            options={'ordering': ('-created_on',)},
        ),
        migrations.AddField(
            model_name='project',
            name='created_on',
            field=models.DateTimeField(auto_now_add=True, db_index=True, null=True),
        ),
        migrations.AddField(
            model_name='project',
            name='description',
            field=models.CharField(blank=True, max_length=255, null=True),
        ),
        migrations.AddField(
            model_name='project',
            name='notes',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='project',
            name='organism',
            field=models.ForeignKey(default=1, on_delete=django.db.models.deletion.CASCADE, to='targets.organism'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='targetcollection',
            name='created_on',
            field=models.DateTimeField(auto_now_add=True, db_index=True, null=True),
        ),
    ]
