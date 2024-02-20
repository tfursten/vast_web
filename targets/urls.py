from django.urls import path

from . import views

app_name = "targets"

urlpatterns = [
    path("", views.index, name="index"),
    path("projects/", views.ProjectListView.as_view(), name='project_list'),
    path("projects/create/", views.ProjectFormView.as_view(), name='project_create'),
    path("projects/<str:pk>/update/", views.ProjectFormView.as_view(), name='project_update'),
    path("projects/<str:pk>/detail/", views.ProjectDetailView.as_view(), name='project_detail'),
    path("projects/<str:pk>/remove/", views.ProjectDeleteView.as_view(), name='project_delete'),
    path("projects/<str:pk>/select-organism/", views.select_organism_list, name="select_organism"),
    path("projects/<str:pk>/add-data/", views.load_project_data, name="add_data"),
    path("projects/<str:pk>/add-loci/", views.select_loci_list, name="select_loci"),
    path("projects/<str:pk>/loci/", views.LocusListView.as_view(), name="locus_list"),
    path("projects/<str:pk>/genomes/", views.GenomeListView.as_view(), name="genome_list"),
    path("projects/<str:pk>/add-genomes/", views.select_genomes_list, name="select_genomes"),
    # path("projects/<str:pk>/loci-data/", views.get_loci_json, name="get_loci_json"),
    # path("projects/<str:pk>/genome-data/", views.get_genomes_json, name="get_genome_json"),
    path("organism/<str:pk>/detail/", views.OrganismDetailView.as_view(), name="organism_detail"),
    path("organisms/", views.OrganismListView.as_view(), name="organism_list"),
    path("organism-data/", views.get_organism_json, name="get_organism_json"),

    


    # path("select-organism/", views.select_new_organism_list, name="select_new_organism"),
    # path("organism/<str:pk>/detail/", views.OrganismDetailView.as_view(), name="organism_detail"),
    # path("organism/<str:organism>/add-loci/", views.add_loci_to_organism, name="add_new_loci"),
    # path("organism/<str:organism>/loci/", views.OrganismLociListView.as_view(), name="organism_loci_list"),
    # path("new-panel/", views.PanelUpdateView.as_view(), name="new_panel"),
    # path("panels/", views.PanelListView.as_view(), name="panel_list"),
    # path("organisms/", views.OrganismListView.as_view(), name="organism_list")

]