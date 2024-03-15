from django.urls import path

from . import views

app_name = "targets"

urlpatterns = [
    path("", views.index, name="index"),
    path("projects/", views.ProjectListView.as_view(), name='project_list'),
    path("projects/create/", views.ProjectFormView.as_view(), name='project_create'),
    path("projects/<str:pk>/update/", views.ProjectUpdateView.as_view(), name='project_update'),
    path("projects/<str:pk>/detail/", views.ProjectDetailView.as_view(), name='project_detail'),
    path("projects/<str:pk>/remove/", views.ProjectDeleteView.as_view(), name='project_delete'),
    path("projects/<str:pk>/add-data/", views.load_allele_data, name="add_data"),
    path("projects/<str:pk>/add-metadata/", views.load_metadata, name="add_metadata"),
    path("projects/<str:pk>/loci/", views.LocusListView.as_view(), name="locus_list"),
    path("projects/<str:pk>/genomes/", views.GenomeListView.as_view(), name="genome_list"),
    path("projects/<str:pk>/metadata-categories/", views.MetadataCategoryListView.as_view(), name="metacat_list"),
    path("projects/<str:pk>/select-loci/", views.select_loci_list, name="select_loci"),
    path("projects/<str:pk>/select-genomes/", views.select_genomes_list, name="select_genomes"),
    path("projects/<str:pk>/select-metadata/", views.select_metacat_list, name="select_metacat"),
    path("projects/<str:pk>/loci-data/", views.get_loci_json, name="get_loci_json"),
    path("projects/<str:pk>/genome-data/", views.get_genomes_json, name="get_genome_json"),
    path("projects/<str:pk>/metadata-data/", views.get_metacat_json, name="get_metacat_json"),
    path("projects/<str:pk>/active-loci-data/<str:selected>/", views.get_loci_json, name="get_loci_json"),
    path("projects/<str:pk>/active-genome-data/<str:selected>/", views.get_genomes_json, name="get_genome_json"),
    path("projects/<str:pk>/active-meta-data/<str:selected>/", views.get_metacat_json, name="get_metacat_json"),
    path("projects/<str:pk>/data/", views.data_list_view, name="data_list"),
    path("projects/<str:pk>/metadata/", views.metadata_list_view, name="metadata_list"),
    # path("projects/<str:pk>/data/get/", views.get_data_json, name="get_data_json"),
    path("projects/<str:pk>/data/remove/", views.data_delete_view, name="data_delete"),
    path("projects/<str:pk>/metadata/remove/", views.metadata_delete_view, name="metadata_delete"),
    path("target-collection/add-target-collection/", views.TargetCollectionFormView.as_view(), name="target_collection_create"),
    path("projects/<str:pk>/add-target-collection/", views.TargetCollectionFormView.as_view(), name="target_collection_create"),
    path("target-collection/<str:pk>/detail/", views.TargetCollectionDetailView.as_view(), name="target_collection_detail"),
    path("target-collection/<str:pk>/update/", views.TargetCollectionUpdateView.as_view(), name="target_collection_update"),
    path("target-collection/list/", views.TargetCollectionListView.as_view(), name="target_collection_list"),
    path("projects/<str:pk>/target-collection/list/", views.TargetCollectionListView.as_view(), name="target_collection_list"),
    path("target-collection/<str:pk>/remove", views.TargetCollectionDeleteView.as_view(), name="target_collection_delete"),
    path("target-collection/<str:pk>/targets/", views.TargetCollectionLocusListView.as_view(), name="target_collection_targets_list"),
    path("target-collection/<str:pk>/targets/edit-targets/", views.select_targets_manual, name="target_collection_add_targets_man"),
    path("target-collection/<str:pk>/targets/add-targets-optimize/", views.TargetCollectionAddTargets.as_view(), name="target_collection_add_targets_opt"),
    path("target-collection/<str:pk>/targets/<int:selected>/", views.get_targets_json, name="get_targets_selected_json"),
    path("target-collection/<str:pk>/targets/all/", views.get_targets_json, name="get_targets_json"),
    path("target-collection/<str:pk>/targets/profile/", views.profile_list_view, name="profile_list"),
    path("target-collection/<str:pk>/resolution/data/", views.get_resolution_figure_html, name="resolution_figure_html"),
    path("target-collection/<str:pk>/resolution/", views.resolution_view, name="resolution_figure"),
    path("target-collection/<str:pk>/tree/data", views.TreeSVGView, name="tree_svg_view"),
    path("target-collection/<str:pk>/tree/", views.tree_view, name="tree_view"),




    


    # path("select-organism/", views.select_new_organism_list, name="select_new_organism"),
    # path("organism/<str:pk>/detail/", views.OrganismDetailView.as_view(), name="organism_detail"),
    # path("organism/<str:organism>/add-loci/", views.add_loci_to_organism, name="add_new_loci"),
    # path("organism/<str:organism>/loci/", views.OrganismLociListView.as_view(), name="organism_loci_list"),
    # path("new-panel/", views.PanelUpdateView.as_view(), name="new_panel"),
    # path("panels/", views.PanelListView.as_view(), name="panel_list"),
    # path("organisms/", views.OrganismListView.as_view(), name="organism_list")

]