import React from "react";

const LoadingOverlay = () => {
    return (
        <div 
            className="loader-overlay" 
            style={{
                position: "fixed", 
                top: 0, 
                left: 0, 
                width: "100%", 
                height: "100%", 
                backgroundColor: "rgba(0, 0, 0, 0.5)", 
                display: "flex", 
                justifyContent: "center", 
                alignItems: "center", 
                zIndex: 1000}}
            >
            <div 
                className="loader" 
                style={{
                    position: "absolute", 
                    top: "50%", 
                    left: "50%", 
                    transform: "translate(-50%, -50%)"
                }} 
            />
        </div>
    );
};

export default LoadingOverlay;