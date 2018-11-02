# Pseudo-random map generation

1. simplex.py

   This class is a simple translation to python of the C++ code from Eliot Eshelman.
   
   It allows for the generation of Perlin noise in 3.5 or 4 dimensions.
   
2. maps.py

   Contains a class to generate pseudo-random maps.
   
   After importing the class and creating an instance, the method help() give you the following guide:

###

 Welcome to my map generator.

 If you are reading this, you already figured out how
 to call this class, so let's skip that part.

 To see a list of the parameters of the map generation
 use the method maps_class.parameter_list().
 The list shows also the methods to change each
 parameter and the valid limits
 The same methods, changing "set" to "get", can be
 used to recall the value of an individual parameter
 The method maps_class.parameter_list() also shows
 the current values.

 To store the map parameters, use
 maps_class.save_parameters(). The default name is
 the parameter "name" plus ".par"
 To load map parameters, use
 maps_class.load_parameters(). The default name is
 the parameter "name" plus ".par"

 To generate a map, use the method
 maps_class.generate_map()

 To generate a wind map, use the method
 maps_class.generate_wind()

 To generate a temperature map, use the method
 maps_class.generate_temperature()

 To generate a moisture map, use the method
 maps_class.generate_moist()

 To generate a biome map, use the method
 maps_class.generate_biome()

 To generate rivers and lakes, use the method
 maps_class.generate_rivers()
 You can use the argument "detailed=True" for a 
 more precise, albeit slower, tracing of rivers.

 Each new element that can be generated requires
 every previous element in the order described above

 You can use maps_class.generate_weather() to generate
 wind, temperature, moisture and biomes once heights
 are done

 You can use maps_class.generate_all() to generate
 every aspect of the map in one go. The mode for the 
 rivers is "defailed=False" with this shortcut

 Everything except the heights requires a full globe map

 The method maps_class.refine() can be used to get a 
 subsection of a map with a much better resolution.
 If the weather and rivers are already calculated, they will
 be also refined, if not, you will not be able to generate them
 afterwards. Be careful, this process is not reversible, 
 it is recommended to save the global map before refining.

 To visualize a map in 2D, use the method
 maps_class.draw(). You can see a list with the
 available projections using the method
 maps_class.help_projection(). To select a
 projection, set the projection parameter to the
 desired short keyword using the set_projection(val)
 method.
 For some projections you may want to rotate the map
 in the longitude axis. Use self.rotate(val) to
 rotate the longitude by val

 To store a generated map, use maps_class.save_map().
 The default name is the parameter "name" plus ".map".
 To load a map, use maps_class.load_map(). The default
 name is the parameter "name" plus ".map".

 To store the map in vtk format, use
 maps_class.save_vtk(). The default name is the
 parameter "name" plus ".vtk"

 Enjoy the maps!
 
# License

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 
###
 
 # TO DO

  - Fix the formatting in help(), did some changes before commiting and did not check the result.
  - Better readme, probably.
  - Improve the visual representation of the wind variable.
  - Revise support for vtk, have not checked since ages ago.

