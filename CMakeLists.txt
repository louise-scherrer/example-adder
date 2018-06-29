CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE(cmake/base.cmake)

SET(PROJECT_NAME example-adder)
SET(PROJECT_DESCRIPTION "This is an example project, to show how to use Gepetto's tools")

SET(PROJECT_URL https://gepgitlab.laas.fr/gepetto/${PROJECT_NAME})

SETUP_PROJECT()

SETUP_PROJECT_FINALIZE()