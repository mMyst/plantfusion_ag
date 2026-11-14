# def simulation(iterations):
#   	### INITIALIZATION ###
#     tagA = "FSPM_a" #<-- these tags are used throughout the script to identify fspm instances
#     tagB = "FSPM_b"
#     indexer = Indexer(global_order=[tagA, tagB], 
#                       FSPMa_names=[tagA],
#                       FSPMb_names=[tagB])

#     planter = Planter()#inputs)

#     fspmA = fspmA_Wrapper(name=tagA) #<-- instances with this tag in 'global_order' will use this wrapper
#     fspmB = fspmB_Wrapper(name=tagB)

#     light = Light_wrapper(#inputs)
#     soil = Soil_wrapper(#inputs)
 
# 	### SIMULATION ###
#     for i in range(1, 1 + iterations):
#       # LIGHT #
#         sceneA = fspmA.light_inputs()
#         sceneB = fspmB.light_inputs()
# 		scenes = indexer.light_scenes_mgmt({tagA : sceneA, tagB : sceneB})
      	
#         light.run(scenes)
      	
#         fspmA.light_results(light)
#         fspmB.light_results(light)

#       # SOIL #
# 		soil_inputsA = fspmA.soil_inputs()
#       	soil_inputsB = fspmB.soil_inputs()
#       	(all_soil_inputs) = indexer.soil_inputs_mgmt({tagA : soil_inputsA, tagB : soil_inputsB})
      
#         soil.run(all_soil_inputs)
           
#         fspmA.soil_results(soil.results)
#       	fspmB.soil_results(soil.results)
     
#       # FSPMs #
#       	fspmA.run()
#         fspmB.run()
      
#     ### END ###
#     fspmA.end()
#     fspmB.end()
#     soil.end()
      
#    if __name__ == "__main__":
#       simulation()