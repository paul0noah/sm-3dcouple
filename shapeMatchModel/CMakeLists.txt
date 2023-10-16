add_subdirectory(constraintsComputation)
add_subdirectory(combinations)
add_subdirectory(energyComputation)

target_sources(${PROJECT_NAME} 
				PRIVATE
					shapeMatchModel.cpp
					shapeMatchModelIlp.cpp
					shapeMatchModel.hpp
	)
