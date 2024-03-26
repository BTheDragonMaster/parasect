import React from "react";

const Navbar = () => {
    return (
        <div
            className="navbar is-fixed-top"
            role="navigation"
            aria-label="main-navigation"
            style={{top: 0, width: "100%", zIndex: 1000, boxShadow: "0 0 10px rgba(0, 0, 0, 0.5)"}}
        >

            {/* Navbar brand. */}
            <div className="navbar-brand">
                <div
                    className="navbar-item" 
                    href="/" 
                    style={{cursor: "pointer"}}
                >
                    {/* Hamburger menu. */}
                    <a
                        role="button"
                        className="navbar-burger"
                        aria-label="menu"
                        aria-expanded="tr"
                        data-target="navbarBasicExample"
                        onClick={() => {
                            const burger = document.querySelector(".navbar-burger");
                            const menu = document.querySelector(".navbar-menu");
                            burger.classList.toggle("is-active");
                            menu.classList.toggle("is-active");
                        }}
                    >
                        <span aria-hidden="true"></span>
                        <span aria-hidden="true"></span>
                        <span aria-hidden="true"></span>
                    </a>
                </div>
            </div>
            
            {/* Navbar menu. */}
            <div id="navbarBasicExample" className="navbar-menu">
                <div className="navbar-start">
                    <a className="navbar-item" href="/">Home</a>
                    <div className="navbar-item has-dropdown is-hoverable">
                        <a className="navbar-link">Apps</a>
                        <div className="navbar-dropdown">
                            <a key={1} className="navbar-item" href="/paras">Paras</a>
                            <a key={2} className="navbar-item" href="/parasect">Parasect</a>
                            <hr className="navbar-divider" />
                            <a 
                                className="navbar-item" 
                                href="https://github.com/BTheDragonMaster/parasect/issues" 
                                target="_blank" 
                                rel="noopener noreferrer"
                            >
                                Report an issue
                            </a>
                        </div>
                    </div> 
                </div>
            </div>
            
        </div>
    );
};

export default Navbar;