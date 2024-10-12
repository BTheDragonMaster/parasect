import React, { useState, useEffect } from "react";

const LoadingOverlay = ({ frame1, frame2 }) => {
    // state to manage which image is currently shown
    const [currentImage, setCurrentImage] = useState(frame1); // assuming the images are named image1.png and image2.png

    // useEffect to toggle the image every second (1000 milliseconds)
    useEffect(() => {
        const intervalId = setInterval(() => {
            setCurrentImage((prevImage) => (prevImage === frame1 ? frame2 : frame1));
        }, 250); // change the interval as needed

        // cleanup function to clear the interval when the component unmounts
        return () => clearInterval(intervalId);
    }, []);

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
            <img 
                width="300px"
                style={{borderRadius: "50%"}}
                src={`${process.env.PUBLIC_URL}/${currentImage}`} 
                alt="Loading..." 
            />
        </div>
    );
};

export default LoadingOverlay;