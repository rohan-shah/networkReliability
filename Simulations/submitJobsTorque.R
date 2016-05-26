source("./generateScenarios.R")
if(!file.exists("./results")) dir.create("results")
if(!exists("maxJobs")) maxJobs <- nrow(scenarios)
submittedJobs <- 0
for(i in 1:nrow(scenarios))
{
	submit <- TRUE
	resultFile <- file.path("results", scenarios[i, "file"])

	shouldSubmit <- FALSE
	if(!file.exists(resultFile)) shouldSubmit <- TRUE
	else
	{
		load(resultFile)
		if(length(results) != scenarios[i, "nReps"]) shouldSubmit <- TRUE
	}
	if(shouldSubmit)
	{
		system2(command = "qsub", args = "submitScriptTorque.sh", env = paste0("SCENARIO_INDEX=", i), wait=TRUE)
		submittedJobs <- submittedJobs + 1
		if(submittedJobs == maxJobs) break
	}
}
