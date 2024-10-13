import React, { useState, useEffect } from 'react';

const Loading = ({ frame1, frame2 }) => {
    // state to manage which image is currently shown
    const [currentImage, setCurrentImage] = useState(frame1);

    // useEffect to toggle the image every 250ms
    useEffect(() => {
        const intervalId = setInterval(() => {
            setCurrentImage((prevImage) => (prevImage === frame1 ? frame2 : frame1));
        }, 500); // adjust the interval timing as needed

        // cleanup function to clear the interval when the component unmounts
        return () => clearInterval(intervalId);
    }, [frame1, frame2]);

    return (
        <div 
            className='loader-container'
            style={{
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                height: '100%', // Take full height of the parent container
                width: '100%'   // Take full width of the parent container
            }}
        >
            <img 
                width='300px'
                src={`${process.env.PUBLIC_URL}/${currentImage}`} 
                alt='Loading...' 
            />
        </div>
    );
};

export default Loading;
