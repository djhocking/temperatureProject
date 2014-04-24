#Index the tile by site location:

indexDaymetTileByLatLon <- function(SiteLat, SiteLon){

  Tile <- ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -74 & SiteLon < -72, 11754,
        ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -72 & SiteLon < -70, 11755,
        ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -70 & SiteLon < -68, 11756,        
        ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -74 & SiteLon < -72, 11934, 
        ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -72 & SiteLon < -70, 11935,
        ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -70 & SiteLon < -68, 11936, #**
        ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -74 & SiteLon < -72, 12114,        
        ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -72 & SiteLon < -70, 12115, 
        ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -70 & SiteLon < -68, 12116, #**      
        ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -68 & SiteLon < -66, 12117, #** 
        ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -72 & SiteLon < -70, 12295, #**       
        ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -70 & SiteLon < -68, 12296, #**      
        ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -68 & SiteLon < -66, 12297, #**    
        print("Tile Error"))))))))))))))

  return(Tile)

}
# ** Assumed, but not explicitly checked.




#Test to see how merging branches works

 A <- 77




#testBranch 2 changes


B <- 3