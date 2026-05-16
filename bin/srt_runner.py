#!python
"""srt_runner.py

Starts the SRT Daemon and/or Dashboard

"""

import argparse
import getpass
import sys
from pathlib import Path
from multiprocessing import Process

from waitress import serve
from srt import config_loader


def run_srt_daemon(configuration_dir, configuration_dict):
    from srt.daemon import daemon as srt_d

    daemon = srt_d.SmallRadioTelescopeDaemon(configuration_dir, configuration_dict)
    daemon.srt_daemon_main()


def run_srt_dashboard(configuration_dir, configuration_dict):
    from srt.dashboard import app as srt_app

    app_server, _ = srt_app.generate_app(configuration_dir, configuration_dict)
    serve(
        app_server,
        host=configuration_dict["DASHBOARD_HOST"],
        port=configuration_dict["DASHBOARD_PORT"],
    )


if __name__ == "__main__":
    # Create the parser
    my_parser = argparse.ArgumentParser(description="Runs the SRT Control Application")

    # Add the arguments
    my_parser.add_argument(
        "--config_dir",
        metavar="config_dir",
        type=str,
        help="The Path to the SRT Config Directory",
        default="~/.srt-config",
    )
    my_parser.add_argument(
        "--config_file_name",
        metavar="config_file_name",
        type=str,
        help="The filename of the Config File to Load",
        default="config.yaml",
    )

    my_parser.add_argument(
        "--dash_only",
        dest="dash_only",
        action="store_true",
        help="Load up the dashboard only",
    )
    my_parser.add_argument(
        "--antenna_user",
        metavar="antenna_user",
        type=str,
        help="Username for antenna control server (W1XMBIGDISH)",
        default=None,
    )
    my_parser.add_argument(
        "--antenna_password",
        metavar="antenna_password",
        type=str,
        help="Password for antenna control server — prefer interactive prompt over this flag",
        default=None,
    )
    # Execute the parse_args() method
    args = my_parser.parse_args()

    # Create Path Objects
    config_dir_path = Path(args.config_dir)
    sky_coords_path = Path(config_dir_path, "sky_coords.csv")
    schema_path = Path(config_dir_path, "schema.yaml")
    config_path = Path(config_dir_path, args.config_file_name)

    dash_only = args.dash_only

    if not config_dir_path.is_dir():
        print("Configuration Directory Does Not Exist")
        print("Please Refer to the Documentation for Creating a Config Directory")
    elif not sky_coords_path.is_file():
        print("Sky Coordinates CSV Not Found")
        print("Please Refer to the Documentation for Creating a sky_coords.csv File")
    elif not schema_path.is_file():
        print("YAML Configuration Schema File Not Found")
        print("Please Put a Copy of the Sample YAML Schema in your Config Directory")
    elif not config_path.is_file():
        print("YAML Configuration File Not Found")
        print("Please Refer to the Documentation for Creating a config.yaml File")
    elif not config_loader.validate_yaml_schema(config_path, schema_path):
        print("YAML Configuration File Invalid")
        print("YAML did not validate against its schema file")
    elif dash_only:
        config_dict = config_loader.load_yaml(config_path)
        dashboard_process = Process(
            target=run_srt_dashboard,
            args=(args.config_dir, config_dict),
            name="SRT-Dashboard",
        )
        dashboard_process.start()
        dashboard_process.join()
    else:
        config_dict = config_loader.load_yaml(config_path)

        if config_dict.get("MOTOR_TYPE") == "W1XMBIGDISH":
            user = args.antenna_user or input("BigDish username: ")
            password = args.antenna_password or getpass.getpass("BigDish password: ")

            from srt.daemon.rotor_control.bigdish_client import BigDishClient
            host = config_dict.get("BIGDISH_HOST", "172.25.15.11")
            port = config_dict.get("BIGDISH_PORT", 1234)
            try:
                probe = BigDishClient(host, port, user, password)
                session_granted = probe.init_session(kick_others=False)
                probe.close()
            except PermissionError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"Could not connect to BigDish server: {e}")
                sys.exit(1)

            kick_others = False
            if not session_granted:
                answer = input("Another session is active. Kick them? [y/N]: ")
                kick_others = answer.strip().lower() == "y"
                if not kick_others:
                    print("Exiting — will not kick existing session.")
                    sys.exit(0)

            config_dict["BIGDISH_USER"] = user
            config_dict["BIGDISH_PASSWORD"] = password
            config_dict["BIGDISH_KICK_OTHERS"] = kick_others

        daemon_process = Process(
            target=run_srt_daemon,
            args=(args.config_dir, config_dict),
            name="SRT-Daemon",
        )
        daemon_process.start()

        if not config_dict["RUN_HEADLESS"]:
            dashboard_process = Process(
                target=run_srt_dashboard,
                args=(args.config_dir, config_dict),
                name="SRT-Dashboard",
            )
            dashboard_process.start()
            dashboard_process.join()
        else:
            daemon_process.join()
