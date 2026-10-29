# Runs CHECKER on CERTIFICATE and fails unless the checker rejects it with its error code (main returns -1),
# i.e. neither accepts it nor crashes
execute_process(COMMAND ${CHECKER} ${CERTIFICATE} RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)

if(NOT result MATCHES "^(255|-1|4294967295)$")
   message(FATAL_ERROR "${CERTIFICATE}: expected a rejection (exit code -1), got ${result}")
endif()
