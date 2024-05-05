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
                className="loader frame1" 
                style={{
                    position: "absolute", 
                    backgroundImage: "url(/paras_loading_frame_1.svg)",
                    // animation: "frame-switch 1s infinite",
                    borderRadius: 0,
                    width: "100%",
                    height: "100%",
                    backgroundSize: "cover",
                    backgroundRepeat: "no-repeat",
                    backgroundPosition: "center"
                }} 
            />
            <div 
                className="loader frame2" 
                style={{
                    position: "absolute", 
                    backgroundImage: "url(/paras_loading_frame_2.svg)",
                    // animation: "frame-switch 1s infinite",
                    borderRadius: 0,
                    width: "100%",
                    height: "100%",
                    backgroundSize: "cover",
                    backgroundRepeat: "no-repeat",
                    backgroundPosition: "center"
                }} 
            />
        </div>
    );
};

export default LoadingOverlay;