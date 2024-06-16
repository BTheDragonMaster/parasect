import React from "react";
import { Link } from "react-router-dom";

const NotFound = () => {

    const logoContainerStyle = {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        flexDirection: 'column', // Stack logo and text vertically
        margin: '20px', // Add some space around each logo container
    };

    return (
        <div 
            style={{
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'center',
                alignItems: 'center',
            }}
        >
            <div style={logoContainerStyle}>
                <div style={{ fontWeight: 'bold', fontSize: '1.5em' }}>
                    404 Not Found
                </div>
                <img 
                    style={{ width: '300px' }}
                    src="/paras_error.png" 
                    alt="PARAS Error"
                />
                <Link to="/">Go back to the home page</Link>
            </div>
        </div>
    );
};

export default NotFound;