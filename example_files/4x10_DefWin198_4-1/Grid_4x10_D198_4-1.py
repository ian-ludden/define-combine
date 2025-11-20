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
        if((x==6) or (x==1 and y==1) or (x==3 and (y==0 or y==3)) or (x==5 and y==2)):
            votes.append(6)
        elif(x>7 or (x==7 and (y==0 or y==3))):
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
sub_dist_cells.append( Polygon([(0, 0), (0, 1), (4,1), (4,0),(0,0)]) )
sub_dist_cells.append( Polygon([(0, 1), (0, 4), (2,4), (2,3), (1,3),(1,1),(0,1)]) )
sub_dist_cells.append( Polygon([(1, 1), (1,2), (3,2), (3,3), (4,3),(4,1),(1,1)]) )
sub_dist_cells.append( Polygon([(1,2), (1,3), (2,3), (2,4), (4,4),(4,3),(3,3),(3,2),(1,2)]) )
sub_dist_cells.append( Polygon([(4, 0), (4,2), (6,2), (6,0),(4,0)]) )
sub_dist_cells.append( Polygon([(4,2), (4,4), (6,4), (6,2), (4,2)]) )
sub_dist_cells.append( Polygon([(6,0), (6,4), (7,4), (7,0), (6,0)]) )
sub_dist_cells.append( Polygon([(7,0), (7,3), (8,3), (8,1), (9,1),(9,0),(7,0)]) )
sub_dist_cells.append( Polygon([(7,3), (7,4), (9,4), (9,1), (8,1),(8,3),(7,3)]) )
sub_dist_cells.append( Polygon([(9,0), (9,4), (10,4), (10,0), (9,0)]) )

# Combine Subdistrict Voting Data - Manual for now
sub_dist_votes = []
sub_dist_votes.append(21)
sub_dist_votes.append(20)
sub_dist_votes.append(21)
sub_dist_votes.append(21)
sub_dist_votes.append(20)
sub_dist_votes.append(21)
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
    # Get representative point
    rep_point = row['geometry'].representative_point()
    # Place vote values
    grid_defined.text(
        rep_point.x, rep_point.y,
        row['votes'],  
        fontsize=10, color=("white" if row['votes']<20 else "black"), ha='center', va='center'
    )