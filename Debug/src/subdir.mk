################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/Faddeeva.cc \
../src/fourier2D.cc 

CPP_SRCS += \
../src/HardCoreBosons.cpp \
../src/PVintegration.cpp \
../src/kernel.cpp 

CC_DEPS += \
./src/Faddeeva.d \
./src/fourier2D.d 

OBJS += \
./src/Faddeeva.o \
./src/HardCoreBosons.o \
./src/PVintegration.o \
./src/fourier2D.o \
./src/kernel.o 

CPP_DEPS += \
./src/HardCoreBosons.d \
./src/PVintegration.d \
./src/kernel.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -O3 -pedantic -Wall -Wextra -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -O3 -pedantic -Wall -Wextra -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


