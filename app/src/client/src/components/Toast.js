import { ToastContainer } from 'react-toastify';

const Toast = () =>  {
    return (
        <ToastContainer
            // set the position of the toast notifications to bottom-right.
            position='bottom-right'
            
            // automatically close the toast after 5000 milliseconds (5 seconds).
            autoClose={5000}
            
            // hide the progress bar for the toast notifications.
            hideProgressBar={true}
            
            // display newer toasts at the top.
            newestOnTop={false}
            
            // disable closing the toast when clicking on it.
            closeOnClick={false}
            
            // set the text direction to left-to-right.
            rtl={false}
            
            // do not pause the toast when focus is lost from the window.
            pauseOnFocusLoss={false}
            
            // disable dragging of toast notifications.
            draggable={false}
            
            // pause the toast when hovering over it.
            pauseOnHover={true}
            
            // customize the toast icons based on the type (success, warning, error, or info).
            icon={({ type }) => {
                if (type === 'success') return '🎉';
                if (type === 'warning') return '⚠️';
                if (type === 'error') return '🚨';
                else return 'ℹ️';
            }}
        />
    );
};

export default Toast;