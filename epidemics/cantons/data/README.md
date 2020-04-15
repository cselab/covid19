# GIS data

<https://shop.swisstopo.admin.ch/en/products/landscape/boundaries3D>

    DATEN/swissBOUNDARIES3D/SHAPEFILE_LV03_LN02/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp

Converted to 2D WGS84 coordinates using QGIS (required by Basemap.readshapefile)

* Layer -> Add Layer -> Add Vector Layer
* Layer -> Save As...
  - CRS: Default CRS ... WGS 84
  - Geometry type: Polygons
  - Include z-dimension: uncheck


# Work-home commuting

<https://www.bfs.admin.ch/bfs/en/home/statistics/mobility-transport/passenger-transport/commuting.assetdetail.8507281.html>

`je-e-11.04.04.05.xlsx`: original

`epidemics/data/files/switzerland_commute_admih_ch.csv`: extraction of

* `Canton initials` for `Commune of residence`
* `Canton initials` for `Commune of work`
* `Number of employed persons`
