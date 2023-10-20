from flask import Flask, render_template
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired


# Create a Flask Instance
app = Flask(__name__)
# TODO: Hide this stuff!!!!
app.config['SECRET_KEY'] = "secretkey"

# Filters
"""
safe - lets you pass html tags
capitalize
lower
upper
title - makes it titlecase
trim - removes trailing spaces
striptags - removes html tags
"""


# Create a Form class

class NamerForm(FlaskForm):
    name = StringField("What's your name?", validators=[DataRequired()])
    submit = SubmitField("Submit")


# Create a route decorator
# A route is a URL - we need to create them. @app refers to created app
@app.route('/')
def index():
    first_name = "Barbara"
    stuff = "This is bold TEXT"
    favourite_pizza = ["gorgonzola", "mushrooms", 41]
    return render_template("index.html",
                           first_name=first_name,
                           stuff=stuff,
                           favourite_pizza=favourite_pizza)


# localhost:5000/user/Barbara
@app.route('/user/<name>')
def user(name):
    return render_template("user.html", name=name)


# Post: move things to the backend
@app.route('/name', methods=["GET", "POST"])
def name():
    name = None
    form = NamerForm()
    if form.validate_on_submit():
        name = form.name.data
        form.name.data = ''
    return render_template("name.html",
                           name=name,
                           form=form)

# Invalid URL

@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404


@app.errorhandler(500)
def internal_server_error(e):
    return render_template("500.html"), 500


