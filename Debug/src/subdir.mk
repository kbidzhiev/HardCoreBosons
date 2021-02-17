################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/HardCoreBosons.cpp \
../src/PVintegration.cpp \
../src/kernel.cpp 

OBJS += \
./src/HardCoreBosons.o \
./src/PVintegration.o \
./src/kernel.o 

CPP_DEPS += \
./src/HardCoreBosons.d \
./src/PVintegration.d \
./src/kernel.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


