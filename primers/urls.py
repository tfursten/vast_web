from django.urls import path

from . import views

app_name = "primers"

urlpatterns = [
    path("primer-collection/", views.PrimerCollectionListView.as_view(), name='primer_collection_list'),
    path("primer-collection/create/", views.PrimerCollectionFormView.as_view(), name='primer_collection_create'),
    path("primer-collection/<str:pk>/update/", views.PrimerCollectionUpdateView.as_view(), name='primer_collection_update'),
    path("primer-collection/<str:pk>/detail/", views.PrimerCollectionDetailView.as_view(), name='primer_collection_detail'),
    path("primer-collection/<str:pk>/remove/", views.PrimerCollectionDeleteView.as_view(), name='primer_collection_delete'),



]