from flask import Flask, Response 

from routes.paras import blueprint_submit_paras
from routes.parasect import blueprint_submit_parasect

app = Flask(__name__)
app.register_blueprint(blueprint_submit_paras)
app.register_blueprint(blueprint_submit_parasect)

@app.errorhandler(404)
def not_found(error) -> Response:
    return app.send_static_file("index.html")

@app.route("/")
def index() -> Response:
    return app.send_static_file("index.html")

def main() -> None:
    app.run(host="localhost", port=4000, debug=True)

if __name__ == "__main__":
    main()