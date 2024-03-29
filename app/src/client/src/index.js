import React, { useState, useEffect } from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route, useLocation } from "react-router-dom";
import { toast } from "react-toastify";

import "react-widgets/scss/styles.scss";
import "./style/main.scss";

import Toast from "./components/Toast";
import Navbar from "./components/Navbar";
import Home from "./components/Home";
import Paras from "./components/Paras";
import Parasect from "./components/Parasect";

function AppRoutes () {
    const location = useLocation();

    const [displayLocation, setDisplayLocation] = useState(location);
    const [transitionStage, setTransitionStage] = useState("fade-in");

    useEffect(() => {
        if (location !== displayLocation) setTransitionStage("fade-out");
    }, [location, displayLocation]);

    return (
        <div
            className={`widget ${transitionStage}`}
            onTransitionEnd={() => {
                if (transitionStage === "fade-out") {
                    setTransitionStage("fade-in");
                    setDisplayLocation(location);
                };
            }}
        >
            <Routes location={displayLocation}>
                <Route path="/" element={<div className="content"><Home /></div>} />
                <Route path="/paras" element={<div className="content"><Paras /></div>} />
                <Route path="/parasect" element={<div className="content"><Parasect /></div>} />
            </Routes>
        </div>
    );
};

function App () {
    return (
        <div>
            <BrowserRouter>
                <Navbar />
                <AppRoutes />
                <Toast />
            </BrowserRouter>
        </div>
    );
};

const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);