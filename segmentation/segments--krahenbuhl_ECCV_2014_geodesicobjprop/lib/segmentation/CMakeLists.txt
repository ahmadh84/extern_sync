add_library( segmentation aggregation.cpp graph.cpp iouset.cpp segmentation.cpp )
target_link_libraries( segmentation util learning imgproc contour )
