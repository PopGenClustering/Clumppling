import logging
import sys
import os

def setup_logger(log_file: str = "logs/app.log"):
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    logger = logging.getLogger()  # root logger
    logger.setLevel(logging.INFO)

    # formatter = logging.Formatter(
    #     '%(asctime)s - %(levelname)s - %(name)s - %(message)s')
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)

    # Stream handler (console)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)

    # Clear previous handlers if re-running in interactive sessions
    if logger.hasHandlers():
        logger.handlers.clear()

    logger.addHandler(fh)
    logger.addHandler(sh)

