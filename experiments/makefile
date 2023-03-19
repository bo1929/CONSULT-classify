MAIN_DIR := $(shell dirname -- "$0";)../
DATASET_PATH := ${MAIN_DIR}/experiments/data
PREDICTIONS_PATH := ${MAIN_DIR}/experiments/predictions
LOGS_PATH := ${MAIN_DIR}/experiments/logs

ensure:
	@mkdir -p ${PREDICTIONS_PATH}
	@mkdir -p ${LOGS_PATH}
	@make consult_classify -C ..

clean:
	@rm -rf ${PREDICTIONS_PATH}
	@rm -rf ${LOGS_PATH}
	@make clean -C ..

sequential:
	@$(MAKE) ensure
	@echo "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "\n===SEQUENTIAL RUN===\n"
	@echo "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@cd .. && \
	./consult_classify \
		-i ./experiments/data \
		-o ./experiments/predictions \
		--taxonomy-path ./experiments/taxonomy/taxonomy-lookup \
		--thread-count 1 > ./experiments/logs/run_numThreads1.out \
	&& cd experiments

parallel:
	@$(MAKE) ensure
	@echo "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "\n===PARALLEL RUN===\n"
	@echo "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@cd .. && \
	read -p "Enter Number of Threads:" numThreads; \
	./consult_classify \
		-i ./experiments/data \
		-o ./experiments/predictions \
		--taxonomy-path ./experiments/taxonomy/taxonomy-lookup \
		--thread-count $${numThreads} > ./experiments/logs/run_numThreads$${numThreads}.out \
	&& cd experiments
