# Coarse country envelopes

`us_boundary.tsv` covers the contiguous United States (48 states).
`uk_boundary.tsv` covers Great Britain, Northern Ireland and northern islands
including Shetland and the Isles of Scilly from the 1:10m source. Remote Rockall
is excluded (retained window: 9 W to 2 E, 49 to 61 N). Check that your final buffered mesh
contains all study locations.
These are small, convex modelling envelopes, not administrative borders or
coastline masks. The UK envelope includes intervening sea and parts of Ireland;
the US envelope excludes Alaska and Hawaii.

Source: Natural Earth Admin 0 Countries, version 5.1.1, public domain.
US: 1:110m; UK: 1:10m. Downloaded 2026-09-10 from
https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip
and https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_countries.zip.
Source documentation:
https://www.naturalearthdata.com/downloads/110m-cultural-vectors/110m-admin-0-countries/.

The bundled coordinate files were prepared locally by constructing a convex
outer envelope from 16 evenly spaced supporting directions in spherical
azimuthal-equidistant coordinates (US origin 38,-97; UK origin 55,-3), then
exporting Lat/Lon degrees. Straight edges deliberately omit coastal detail.
The analysis example projects these coordinates into the model's frame and
buffers the convex envelope by a fixed 50 km (US) or 5 km (UK), independently
of the nominal correlation range.
