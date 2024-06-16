import React, {useState} from "react";

const Home = () => {
    const [isHoveredParas, setIsHoveredParas] = useState(false);
    const [isHoveredParasect, setIsHoveredParasect] = useState(false);

    const logoContainerStyle = {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        flexDirection: 'column', // Stack logo and text vertically
        margin: '20px', // Add some space around each logo container
    };

    const logoStyleParas = {
        width: '300px',
        transform: isHoveredParas ? 'scale(1.1)' : 'none', 
        transition: 'transform 0.2s',
        cursor: 'pointer',
        borderRadius: '50%',
        marginBottom: '10px',
    };

    const logoStyleParasect = {
        width: '300px',
        transform: isHoveredParasect ? 'scale(1.1)' : 'none',
        transition: 'transform 0.2s',
        cursor: 'pointer',
        borderRadius: '50%',
        marginBottom: '10px',
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
            <div style={{ fontSize: '1.5em', marginBottom: '20px' }}>
                Welcome to Paras! Pick an app from below or the navigation menu above to get started.
            </div>
            <div style={{ fontSize: '1.5em', marginBottom: '20px' }}>
                Read more about PARAS and PARASECT in our <a href="/publication">publication</a>.
            </div>
            <div 
                style={{
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    flexDirection: 'row', // Align logos horizontally
                }}
            >
                <div style={logoContainerStyle}>
                    <img 
                        src="/paras.png" 
                        alt="PARAS" 
                        style={logoStyleParas} 
                        onMouseEnter={() => setIsHoveredParas(true)}
                        onMouseLeave={() => setIsHoveredParas(false)}
                        onClick={() => window.location.href = '/paras'}
                    />
                    <div style={{ fontWeight: 'bold', fontSize: '1.5em' }}>
                        PARAS
                    </div>
                </div>
                <div style={logoContainerStyle}>
                    <img 
                        src="/parasect.png" 
                        alt="PARASECT" 
                        style={logoStyleParasect} 
                        onMouseEnter={() => setIsHoveredParasect(true)}
                        onMouseLeave={() => setIsHoveredParasect(false)}
                        onClick={() => window.location.href = '/parasect'}
                    />
                    <div style={{ fontWeight: 'bold', fontSize: '1.5em' }}>              
                        PARASECT
                    </div>
                </div>
            </div>
        </div>
    );
};

export default Home;