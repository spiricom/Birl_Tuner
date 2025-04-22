cmake_minimum_required(VERSION 3.25)

# This tells cmake we have goodies in the /cmake folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Modern concise way to add dependencies to your project

# Configures universal binaries and decides which version of macOS to support

# Couple tweaks that IMO should be JUCE defaults

# Change me!
# This is the internal name of the project and the name of JUCE's shared code target
# Note: This cannot have spaces (it may be 2024, but you can't have it all!)
# Worry not, JUCE's PRODUCT_NAME can have spaces (and is what DAWs display)
set(PROJECT_NAME "Pamplejuce")

# Worry not, JUCE's PRODUCT_NAME can have spaces (and is what DAWs will display)
# You can also just have it be the same thing as PROJECT_NAME
# You may want to append the major version on the end of this (and PROJECT_NAME) ala:
#   set(PROJECT_NAME "MyPlugin_v${MAJOR_VERSION}")
# Doing so enables major versions to show up in IDEs and DAWs as separate plugins
# allowing you to change parameters and behavior without breaking existing user projects
set(PRODUCT_NAME "Pamplejuce Demo")

# Change me! Used for the MacOS bundle name and Installers
set(COMPANY_NAME "Pamplejuce Company")

# Change me! Used for the MacOS bundle identifier (and signing)
set(BUNDLE_ID "com.pamplejuce.pamplejuce")

# Change me! Set the plugin formats you want built
# Valid choices: AAX Unity VST VST3 AU AUv3 Standalone
set(FORMATS Standalone )

# For simplicity, the name of the CMake project is also the name of the target
project(${PROJECT_NAME} VERSION ${CURRENT_VERSION})


# JUCE is setup as a submodule in the /JUCE folder
# Locally, you must run `git submodule update --init --recursive` once
# and later `git submodule update --remote --merge` to keep it up to date
# On Github Actions, this is done as a part of actions/checkout
add_subdirectory(JUCE)
add_subdirectory(third_party/LEAF)
include_directories(third_party/LEAF/leaf)

# 1) Collect just the DSP/tuning code you need
file(GLOB_RECURSE CORE_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/*.h"
)

# 2) Build a static lib from it
add_library(BirlCore STATIC ${CORE_SOURCES})
target_link_libraries(BirlCore
        PUBLIC
        LEAF
        juce::juce_graphics
        juce::juce_core
        juce::juce_dsp
        juce::juce_audio_basics

)

# 3) Your sim executable
add_executable(BirlSimulation
        source/TuningOpt.cpp
)
target_link_libraries(BirlSimulation
        PRIVATE BirlCore
)
target_include_directories(BirlSimulation
        PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/source"
)

# 4) (Optional) always drop the binary in build/bin/
set_target_properties(BirlSimulation PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
# Add CLAP format




