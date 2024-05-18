import os
import datetime

# Absolute path to the logs folder in the project
PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "logs"))

class Logger:
    def __init__(self, path: str = PATH) -> None:
        self.path = path
        self.date = datetime.datetime
        self.log("Log opened", "OPEN LOG")

    def log(self, message: str = "",report_type: str = "DEBUG") -> None:
        with open(os.path.join(self.path, "vgf-scattering.log"), "a") as log_file:
            now = self.date.now().strftime("%d-%b-%Y: %H:%M:%S %z")
            log_file.write(f"{now}[{report_type}]: {message}\n")
