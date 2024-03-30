import React from "react";

const Modal = ({ children, closeModal, modalState, title }) => {
    
    if(!modalState) {
      return null;
    };
    
    return(
        <div className="modal is-active">
            <div 
                className="modal-background" 
                onClick={closeModal} 
            />
                <div 
                    className="modal-card" 
                    style={{
                        borderRadius: "10px", 
                        width: "90%", 
                        height: "90%", 
                        paddingTop: "50px"
                    }}
                >
                    <header className="modal-card-head">
                        <p className="modal-card-title">
                            {title}
                        </p>
                        <button 
                            className="delete" 
                            onClick={closeModal} 
                        />
                    </header>
                <section className="modal-card-body">
                    <div className="content">
                        {children}
                    </div>
                </section>
            </div>
        </div>
    );
};

export default Modal;