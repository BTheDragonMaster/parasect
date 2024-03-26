import { ToastContainer } from "react-toastify";

const Toast = () =>  {
    return (
        <ToastContainer
            // Set the position of the toast notifications to bottom-right.
            position="bottom-right"
            
            // Automatically close the toast after 5000 milliseconds (5 seconds).
            autoClose={5000}
            
            // Hide the progress bar for the toast notifications.
            hideProgressBar={true}
            
            // Display newer toasts at the top.
            newestOnTop={false}
            
            // Disable closing the toast when clicking on it.
            closeOnClick={false}
            
            // Set the text direction to left-to-right.
            rtl={false}
            
            // Do not pause the toast when focus is lost from the window.
            pauseOnFocusLoss={false}
            
            // Disable dragging of toast notifications.
            draggable={false}
            
            // Pause the toast when hovering over it.
            pauseOnHover={true}
            
            // Customize the toast icons based on the type (success, warning, error, or info).
            icon={({ type }) => {
                if (type === "success") return "🎉";
                if (type === "warning") return "⚠️";
                if (type === "error") return "🚨";
                else return "ℹ️";
            }}
        />
    );
};

export default Toast;