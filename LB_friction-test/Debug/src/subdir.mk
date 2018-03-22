################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/LB_Analyze.c \
../src/LB_GUI.c \
../src/LB_Initialize.c \
../src/LB_Simulation.c \
../src/LB_friction-test.c 

OBJS += \
./src/LB_Analyze.o \
./src/LB_GUI.o \
./src/LB_Initialize.o \
./src/LB_Simulation.o \
./src/LB_friction-test.o 

C_DEPS += \
./src/LB_Analyze.d \
./src/LB_GUI.d \
./src/LB_Initialize.d \
./src/LB_Simulation.d \
./src/LB_friction-test.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -std=gnu99 -I"/home/clark/school/Lattice Boltzmann/c/newgraph/include" -O0 -g -Wall -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


