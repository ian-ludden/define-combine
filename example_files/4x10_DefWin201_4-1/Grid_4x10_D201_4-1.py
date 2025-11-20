import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

# Define the extent of your grid (xmin, ymin, xmax, ymax)
xmin, ymin, xmax, ymax = 0, 0, 10, 4

# Define the cell size
cell_size = 1

# Create lists to store grid cell polygons
grid_cells = []

# Generate grid cells
for x in np.arange(xmin, xmax, cell_size):
    for y in np.arange(ymin, ymax, cell_size):
        polygon = Polygon([(x, y), (x + cell_size, y), (x + cell_size, y + cell_size), (x, y + cell_size)])
        grid_cells.append(polygon)

# Assign Vote Counts
votes = []
for x in range(0, 10):
    for y in range(0, 4):
        if((x==0 and y==0) or (x==2 and (y==1 or y==2)) or (x==4 and y==3) or (x>3 and x<7 and y==2) or (x>4 and x<8 and y==1) or (x==7 and y==0)):
            votes.append(6)
        elif((x>4 and y==3) or (x>7 and not (x==9 and y==0))):
            votes.append(4)
        else: votes.append(5)
    
# Define data for file creation
data = {"votes": votes, "geometry": grid_cells}

# Create a GeoDataFrame from the grid cells
grid_gdf = gpd.GeoDataFrame(data, crs="EPSG:4326") # Replace EPSG:4326 with your desired CRS

# Save the GeoDataFrame as a shapefile
grid_gdf.to_file("basic_grid.shp")

# Plot 
grid_predefine = grid_gdf.plot(legend=True, column="votes", edgecolor="black", cmap="bwr")
grid_predefine.set_axis_off()

for idx, row in grid_gdf.iterrows():
    # Get representative point
    rep_point = row['geometry'].representative_point()
    # Place vote values
    grid_predefine.text(
        rep_point.x, rep_point.y,
        row['votes'],
        fontsize=10, color=("white" if row['votes']<5 else "black"), ha='center', va='center'
    )

# Create Subdistrict Grid - Manully generated for now
sub_dist_cells = []
sub_dist_cells.append( Polygon([(0, 0), (0, 1), (1,1), (1,2), (2,2),(2,1),(3,1),(3,0),(0,0)]) )
sub_dist_cells.append( Polygon([(0, 1), (0, 4), (2,4), (2,3), (1,3),(1,1),(0,1)]) )
sub_dist_cells.append( Polygon([(1, 2), (1,3), (2,3), (2,4), (4,4),(4,3),(3,3),(3,2),(1,2)]) )
sub_dist_cells.append( Polygon([(2,1), (2,2), (3,2), (3,3), (4,3),(4,2),(5,2),(5,1),(2,1)]) )
sub_dist_cells.append( Polygon([(3, 0), (3, 1), (7,1), (7,0)]) )
sub_dist_cells.append( Polygon([(4,2), (4,4), (5,4), (5,3), (7,3),(7,2),(4,2)]) )
sub_dist_cells.append( Polygon([(5,1), (5,2), (8,2), (8,0), (7,0),(7,1),(5,1)]) )
sub_dist_cells.append( Polygon([(5,3), (5,4), (8,4), (8,2), (7,2),(7,3),(5,3)]) )
sub_dist_cells.append( Polygon([(8, 0), (8,3), (9,3), (9,1), (10,1),(10,0),(8,0)]) )
sub_dist_cells.append( Polygon([(8,3), (8,4), (10,4), (10,1), (9,1),(9,3),(8,3)]) )

# Combine Subdistrict Voting Data - Manual for now
sub_dist_votes = []
sub_dist_votes.append(21)
sub_dist_votes.append(20)
sub_dist_votes.append(21)
sub_dist_votes.append(21)
sub_dist_votes.append(20)
sub_dist_votes.append(24)
sub_dist_votes.append(24)
sub_dist_votes.append(17)
sub_dist_votes.append(17)
sub_dist_votes.append(16)

# Create file for subdistrict graph data
sub_dist_data = {"votes": sub_dist_votes, "geometry": sub_dist_cells}
sub_dist_gdf = gpd.GeoDataFrame(sub_dist_data, crs="EPSG:4326") # Replace EPSG:4326 with your desired CRS
sub_dist_gdf.to_file("subdistrict_grid.shp")

# Plot - TODO: remove axes, add vote value labels at representative points
grid_defined = sub_dist_gdf.plot(legend=True, column="votes", edgecolor="black", cmap="bwr")
grid_defined.set_axis_off()

for idx, row in sub_dist_gdf.iterrows():
    # Get centroid coordinates
    centroid = row['geometry'].representative_point()
    # Place text (e.g., country name)
    grid_defined.text(
        centroid.x, centroid.y,
        row['votes'],  # Column to label
        fontsize=10, color=("white" if row['votes']<20 else "black"), ha='center', va='center'
    )