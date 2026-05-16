"""system_page.py

Function for Generating System Page and Creating Callback

"""

try:
    from dash import dcc
except:
    import dash_core_components as dcc

try:
    from dash import html
except:
    import dash_html_components as html

import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from urllib.parse import quote as urlquote
from datetime import datetime
from pathlib import Path


def generate_layout(config=None):
    """Generates the Basic Layout for the System Page

    Returns
    -------
    System Page Layout
    """
    is_bigdish = config is not None and config.get("MOTOR_TYPE") == "W1XMBIGDISH"

    bigdish_panel = html.Div(
        [
            html.Div(
                [
                    html.H4(
                        "Antenna Session",
                        style={"text-align": "center"},
                    ),
                    dcc.Markdown(id="bigdish-status"),
                    html.Div(
                        [
                            dbc.Button(
                                "Connect",
                                id="bigdish-connect-btn",
                                color="primary",
                                className="me-2",
                                style={"margin-right": "8px"},
                            ),
                            dbc.Button(
                                "Kick & Connect",
                                id="bigdish-kick-btn",
                                color="danger",
                            ),
                        ],
                        style={"text-align": "center", "margin-top": "8px"},
                    ),
                    html.Div(id="bigdish-connect-result"),
                ],
                className="pretty_container four columns",
            ),
        ],
        className="flex-display",
        style={"justify-content": "center", "margin": "5px"},
    ) if is_bigdish else html.Div()

    layout = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [],
                        className="one-third column",
                    ),
                    html.Div(
                        [
                            html.H4(
                                "SRT System Page",
                                style={"margin-bottom": "0px", "text-align": "center"},
                            ),
                        ],
                        className="one-third column",
                        id="title",
                    ),
                    html.Div([], className="one-third column", id="button"),
                ],
                id="header",
                className="row flex-display",
                style={"margin-bottom": "25px"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.H4(
                                "Emergency Contact Info",
                                id="text-contact",
                                style={"text-align": "center"},
                            ),
                            dcc.Markdown(id="emergency-contact-info"),
                        ],
                        className="pretty_container four columns",
                    ),
                    html.Div(
                        [
                            html.H4(
                                id="text-queue-status", style={"text-align": "center"}
                            ),
                            dcc.Markdown(id="command-display"),
                        ],
                        className="pretty_container four columns",
                    ),
                    html.Div(
                        [
                            html.H4(
                                "Recordings",
                                id="text-recordings",
                                style={"text-align": "center"},
                            ),
                            html.Div(
                                id="recordings-list",
                                style={
                                    "height": 150,
                                    "overflow": "hidden",
                                    "overflow-y": "scroll",
                                },
                            ),
                        ],
                        className="pretty_container four columns",
                    ),
                ],
                className="flex-display",
                style={"justify-content": "center"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.H4(
                                "Message Logs",
                                style={"text-align": "center"},
                            ),
                            html.Div(
                                id="message-logs",
                                style={
                                    "height": 200,
                                    "overflow": "hidden",
                                    "overflow-y": "scroll",
                                },
                            ),
                        ],
                        className="pretty_container twelve columns",
                    ),
                ],
                className="flex-display",
                style={"justify-content": "center", "margin": "5px"},
            ),
            bigdish_panel,
        ]
    )
    return layout


def register_callbacks(app, config, status_thread, command_thread=None):
    """Registers the Callbacks for the System Page

    Parameters
    ----------
    app : Dash Object
        Dash Object to Set Up Callbacks to
    config : dict
        Contains All Settings for Dashboard / Daemon
    status_thread : Thread
        Thread for Getting Status from Daemon
    command_thread : Thread, optional
        Thread for Sending Commands to Daemon

    Returns
    -------
    None
    """
    is_bigdish = config.get("MOTOR_TYPE") == "W1XMBIGDISH"

    @app.callback(
        Output("emergency-contact-info", "children"),
        [Input("interval-component", "n_intervals")],
    )
    def update_contact_info(n):
        status = status_thread.get_status()
        if status is None or "emergency_contact" not in status:
            return ""
        status = status["emergency_contact"]
        status_string = f"""
         - Name: {status["name"]}
         - Email: {status["email"]}
         - Phone Number: {status["phone_number"]}
        """
        return status_string

    @app.callback(
        Output("message-logs", "children"),
        [Input("interval-component", "n_intervals")],
    )
    def update_message_logs(n):
        status = status_thread.get_status()
        if status is None or "error_logs" not in status:
            return ""
        status = status["error_logs"]
        children = [
            html.P(
                f"{datetime.fromtimestamp(log_time).strftime('%Y-%m-%d %H:%M:%S')}: {log_txt}"
            )
            for log_time, log_txt in status
        ]
        return html.Div(children=children)

    @app.callback(
        Output("text-queue-status", "children"),
        [Input("interval-component", "n_intervals")],
    )
    def update_command_queue_display(n):
        status = status_thread.get_status()
        if status is None or (
            status["queue_size"] == 0 and status["queued_item"] == "None"
        ):
            return "SRT Inactive"
        return "SRT in Use!"

    @app.callback(
        Output("command-display", "children"),
        [Input("interval-component", "n_intervals")],
    )
    def update_command_display(n):
        status = status_thread.get_status()
        if status is None:
            return ""
        current_cmd = status["queued_item"]
        queue_size = status["queue_size"]
        status_string = f"""
        ##### Command Queue Status
         - Running Command: {current_cmd}
         - {queue_size} More Commands Waiting in the Queue
        """
        return status_string

    @app.callback(
        Output("recordings-list", "children"),
        [Input("interval-component", "n_intervals")],
    )
    def update_output(n):
        """Save uploaded files and regenerate the file list."""

        files = [
            file.name
            for file in Path(config["SAVE_DIRECTORY"]).expanduser().glob("*")
            if file.is_file()
        ]
        folders = [
            file.name
            for file in Path(config["SAVE_DIRECTORY"]).expanduser().glob("*")
            if file.is_dir()
        ]
        if len(files) == 0:
            return [html.Li("No files yet!")]
        else:
            if config["DASHBOARD_DOWNLOADS"]:
                return [
                    html.Li(html.A(filename, href=f"/download/{urlquote(filename)}"))
                    for filename in files
                ] + [html.Li(html.A(foldername)) for foldername in folders]
            else:
                return [html.Li(html.A(filename)) for filename in (files + folders)]

    if is_bigdish:
        @app.callback(
            Output("bigdish-status", "children"),
            [Input("interval-component", "n_intervals")],
        )
        def update_bigdish_status(n):
            status = status_thread.get_status()
            if status is None:
                return "Daemon not connected"
            busy = status.get("bigdish_session_busy", False)
            if busy:
                return "**Session busy** — another client is active"
            return "**Session active**"

        @app.callback(
            Output("bigdish-connect-result", "children"),
            [
                Input("bigdish-connect-btn", "n_clicks"),
                Input("bigdish-kick-btn", "n_clicks"),
            ],
            prevent_initial_call=True,
        )
        def handle_bigdish_buttons(connect_clicks, kick_clicks):
            from dash import ctx
            if command_thread is None:
                return "No command thread available"
            triggered = ctx.triggered_id
            if triggered == "bigdish-connect-btn":
                command_thread.add_to_queue("bigdish_connect")
            elif triggered == "bigdish-kick-btn":
                command_thread.add_to_queue("bigdish_connect kick")
            return ""
