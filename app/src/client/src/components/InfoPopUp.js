import React, { useState } from "react";
import { MdInfoOutline } from "react-icons/md";

const InfoPopUp = ({ infoText }) => {
    const [isVisible, setIsVisible] = useState(false);

    if (!infoText) {
        infoText = "No information available.";
    };
  
    return (
        <div style={{display: "flex", alignItems: "center", position: "relative"}}>
            <MdInfoOutline
                onMouseEnter={() => setIsVisible(true)}
                onMouseLeave={() => setIsVisible(false)}
                style={{cursor: "pointer", marginLeft: "5px"}}
            />
            <div 
                style={{
                    visibility: isVisible ? "visible" : "hidden",
                    position: "absolute",
                    marginLeft: "25px",
                    padding: "15px",
                    backgroundColor: "lightgray",
                    borderRadius: "5px",
                    boxShadow: "0 0 5px rgba(0,0,0,0.3)",
                    zIndex: 1,
                    width: "210px",
                    wordWrap: "break-word"
                }}
            >
                {infoText}
            </div>
        </div>
    );
};
  
export default InfoPopUp;