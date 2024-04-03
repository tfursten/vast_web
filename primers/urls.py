from django.urls import path

from . import views

app_name = "primers"

urlpatterns = [
    path("primer-collections/", views.PrimerCollectionListView.as_view(), name='primer_collection_list'),
    path("primer-collections/create/", views.PrimerCollectionFormView.as_view(), name='primer_collection_create'),
    path("primer-collections/<str:pk>/update/", views.PrimerCollectionUpdateView.as_view(), name='primer_collection_update'),
    path("primer-collections/<str:pk>/detail/", views.PrimerCollectionDetailView.as_view(), name='primer_collection_detail'),
    path("primer-collections/<str:pk>/remove/", views.PrimerCollectionDeleteView.as_view(), name='primer_collection_delete'),
    path("amplicon-targets/", views.AmpliconTargetListView.as_view(), name='amplicon_target_list'),
    path("amplicon-targets/create/", views.AmpliconTargetFormView.as_view(), name='amplicon_target_create'),
    path("amplicon-targets/<str:pk>/update/", views.AmpliconTargetUpdateView.as_view(), name='amplicon_target_update'),
    path("amplicon-targets/<str:pk>/detail/", views.AmpliconTargetDetailView.as_view(), name='amplicon_target_detail'),
    path("amplicon-targets/<str:pk>/remove/", views.AmpliconTargetDeleteView.as_view(), name='amplicon_target_delete'),


]