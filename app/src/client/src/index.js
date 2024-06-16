import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route } from "react-router-dom";

import "react-widgets/scss/styles.scss";
import "./style/main.scss";

import Toast from "./components/Toast";
import Navbar from "./components/Navbar";
import Home from "./components/Home";
import Paras from "./components/Paras";
import Parasect from "./components/Parasect";
import NotFound from "./pages/NotFound";

function AppRoutes () {
    return (
        <div>
            <Routes>
                <Route path="/" element={<div style={{padding: "30px"}}><Home /></div>} />
                <Route path="/paras" element={<div style={{padding: "20px"}}><Paras /></div>} />
                <Route path="/parasect" element={<div style={{padding: "20px"}}><Parasect /></div>} />
                <Route path="*" element={<div style={{padding: "30px"}}><NotFound /></div>} />
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